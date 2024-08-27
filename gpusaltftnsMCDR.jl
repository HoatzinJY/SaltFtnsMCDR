#import relevant packages and functions 
using Pkg
using Oceananigans
using OceanBioME
using TimesDates
using CairoMakie
using Printf
using Oceananigans, SeawaterPolynomials.TEOS10
using Oceananigans.Models: seawater_density
using SeawaterPolynomials
using GibbsSeaWater
using NetCDF
using DSP
using MPIPreferences
using Oceanostics
using Oceanostics: _νᶜᶜᶜ
using Oceananigans.Operators
using Oceananigans.TurbulenceClosures: viscous_flux_wx
using Oceananigans: NonhydrostaticModel, BuoyancyField
using Oceananigans.BuoyancyModels: buoyancy
using Adapt
using Oceananigans.Grids: NegativeZDirection, validate_unit_vector
using Oceananigans.BuoyancyModels: buoyancy_perturbationᶜᶜᶜ, get_temperature_and_salinity, ρ′


#constants for ease of code use 
const day = 86400;
const hour = 3600;
const minute = 60;

#search #SIMILITUDE to change things for similutde theories 
#assume base SI units for anything that is unspecified (including time, in s)
#base/original parameters used to develop model are commented on the side 



"""NAME"""
#this section creates folder that saves simulation data to. it assumes that you have created a 
#folder in the current directory titled "Trials" to store the data in 
trial_name = "SIMILITUDE to big longer"
mkdir(joinpath("Trials", (trial_name)))
pathname = joinpath("Trials", (trial_name))
filename = joinpath(pathname, "data")



"""COMPUTER parameters"""
#to ensure you don't overload on the grid points 
const GPU_memory = 12



"""SIMULATION RUN INFORMATION"""
simulation_duration = 1hour 
run_duration = 12hour #wall time
output_interval = 10 



"""DOMAIN SIZE & SETUP"""
const domain_x = 0.7;#SIMILITUDE #orig 0.7
const domain_z = 2; #SIMILITUDE #orig 2
const x_center = domain_x/2;
const max_grid_spacing = 0.001 #SIMILITUDE orig 0.001. Somewhat of a misnomer, this will effectively be the grid spacing



"""PIPE PROPERTY INFORMATION"""
#to store physical properties of pipe wall materials 
struct PipeWallData 
    thermal_diffusivity :: Float64
    thickness :: Float64
end
const pipe_radius = 0.02 #SIMILITUDE #orig 0.02
const pipe_length = 1 #SIMILITUDE #oric 1
const pipe_top_depth = 0.5 #SIMILITUDE #orig 0.5 
const pipe_wall_thickness_intended = 0.005 #Similitude, 0.005 orig
const pipe_bottom_depth = pipe_top_depth + pipe_length
const wall_material_ρ = 8900
const wall_material_cₚ = 376.812 
const wall_material_k = 50.208
const pipe_data = PipeWallData(wall_material_k/(wall_material_cₚ*wall_material_ρ), pipe_wall_thickness_intended)



"""CONSTANTS AND PHYSICAL PROPERTIES"""
#utility constants
const g = 9.806;
const geopotential_height = 0;
const CFL = 0.2;
#diffusion information
struct SeawaterDiffusivityData
    ν :: Float64
    T :: Float64
    S :: Float64
end
const eddy_horizontal_diffusivity = 5e2 #not used
const eddy_vertical_diffusivity = 1e-5 #not used
const sw_viscosity_molecular =  1.05e-6 #SIMILITUDE
const sw_T_diffusivity_molecular = 1.46e-7 #SIMILITUDE, note that this is 1e-5 per experimental data Zhang 2004
const sw_S_diffusivity_molecular = 1.3e-9 
const sw_diffusivity_data = SeawaterDiffusivityData(sw_viscosity_molecular, sw_T_diffusivity_molecular, sw_S_diffusivity_molecular)



"""SET BACKGROUND TEMPERATURE AND SALINITY GRADIENTS"""
const T_top = 21.67;
const T_bot = 11.86;
const S_bot = 34.18;
const S_top = 35.22;
const delta_z = 200; #distance between top & bottom measurements
const delta_z_multiplier = 1; #SIMILITUDE, smaller => steeper gradeint 

#=
TWaterColumn & SWaterColumn give you the background water column temperature & salinity properties 
respectively at depth z. Note that it follows model convention for z, such that it expects 
surface to have depth 0 & increasingly negative depths below. Currently, it creates a linear gradient
from the information provided above for the length of the pipe, but sets uniform properties for the region 
above and below the pipe (the same properties as the very top/bottom of the pipe). 

These functions can be edited to interpolate the information from a netcdf file. 
=#
function TWaterColumn(z)
    if (z > -pipe_top_depth)
        return T_top - ((T_bot - T_top) / delta_z)*(-pipe_top_depth)
    elseif (z < -pipe_bottom_depth)
        return T_top - ((T_bot - T_top) / delta_z)*(-pipe_bottom_depth)
    else
        return T_top - ((T_bot - T_top) / delta_z)z
    end
end

function SWaterColumn(z)
    if (z > -pipe_top_depth)
        return S_top - ((S_bot - S_top) / (delta_z * delta_z_multiplier))*(-pipe_top_depth)
    elseif (z < -pipe_bottom_depth)
        return S_top - ((S_bot - S_top) / (delta_z * delta_z_multiplier))*(-pipe_bottom_depth)
    else
        return S_top - ((S_bot - S_top) / (delta_z * delta_z_multiplier))z
    end
end

TWaterColumn(x, z) = TWaterColumn(z)
TWaterColumn(x, y, z) = TWaterColumn(z)
SWaterColumn(x, z) = SWaterColumn(z)
SWaterColumn(x, y, z) = SWaterColumn(z)
TWaterColumnBCS(y, z, t) = TWaterColumn(z)
TWaterColumnBCS(z, t) = TWaterColumn(z)
SWaterColumnBCS(y, z, t) = SWaterColumn(z)
SWaterColumnBCS(z, t) = SWaterColumn(z)


"""COMPUTATIONAL FUNCTIONS"""
#roundup takes in two numbers num and base, and then rounds num up to the 
#nearest multiple of base 
roundUp(num::Float64, base) = ceil(Int, num / base) * base

#findNearest finds the nearest value in array A to given value x. It returns this 
#information as a tuple (minimum difference, [cartesian index of found closest match])
findNearest(A::AbstractArray, x) = findmin(abs(A-x))

#getXYZ returns the [x, y, z] spatial coordinates for given discrete indices of a grid. 
#the function takes indices i, j, k, a grid, and the locations (face, center, etc) the grid
#is evaluated at 
function getXYZ(i, j, k, grid, locs)
    myNodes = nodes(grid, locs, reshape = true)
    return[myNodes[1][i], myNodes[2][j], myNodes[3][k]]
end

#getMaxandMin returns a tuple (minimum val, maximum val) of the values in a certain 
#field time series for all times in the series. The function takes in numPoints (number of times
#in time series) & dataSeries (the time series)
function getMaxAndMin(numPoints, dataSeries)
    myMax = maximum(interior(dataSeries[1], :, 1, :))
    myMin = minimum(interior(dataSeries[1], :, 1, :))
    for i in 2:numPoints
        myMax = max(maximum(interior(dataSeries[i], :, 1, :)),myMax )
        myMin = min(minimum(interior(dataSeries[i], :, 1, :)), myMin)
    end
    return (myMin, myMax)
end

#getMaxAndMinArr returns a tuple (minimum val, maximum val) for the array that is passed in 
function getMaxAndMinArr(arr)
    myMax = maximum(arr)
    myMin = minimum(arr)
    return (myMin, myMax)
end

#getMaxAndMinEnd returns a tuple (minimum val, maximum val) for the last dataset that is saved in 
#the field time series dataseries that is passed in. This is generally useful for plotting, as 
#often the largest values are seen at the end
function getMaxAndMinEnd(dataSeries)
    myMax = maximum(interior(dataSeries[end], : , 1, :))
    myMin = minimum(interior(dataSeries[end], : , 1, :))
    return (myMin, myMax)
end

#getMaskedAverage returns the average of a field over a masked region. The function multiplies 
#each point in the field with the value of the mask evaluated at that point. The function takes in 
#mask(x, z) function, a field, and locs on the field that it is evaluated at. The function could be 
#in theory used to get weighted averages, if the mask acts as the weight. 
function getMaskedAverage(mask::Function, field, locs)
    fieldNodes = nodes(field.grid, locs, reshape = true)
    sum = 0;
    count = 0;
    for i in eachindex(fieldNodes[1])
        for j in eachindex(fieldNodes[3])
            sum += mask(fieldNodes[1][i], fieldNodes[3][j]) * field.data[i, 1, j]
            count += 1*mask(fieldNodes[1][i], fieldNodes[3][j])
        end
    end
    return sum/count
end
getMaskedAverageTracer(mask::Function, field) = getMaskedAverage(mask::Function, field, (Center(), Center(), Center()))
getMaskedAverageVelocity(mask::Function, field) = getMaskedAverage(mask::Function, field, (Face(), Face(), Face()))

#getMaskedAverageAtZ does the same thing as above except it only evaluates the average of 
#the masked area at one z location. The function takes in a mask(x, z), a field with data, locs 
#that the data is saved at, the spatial zCoord (in m), and the spacing in the z. It will take the average 
#of the two grid lines that "sandwich" the given zCoord
function getMaskedAverageAtZ(mask::Function, field, locs, zCoord, z_spacings)
    fieldNodes = nodes(field.grid, locs, reshape = true)
    sum = 0;
    count = 0;
    for i in eachindex(fieldNodes[1])
        for j in eachindex(fieldNodes[3])
            if (abs(fieldNodes[3][j] - zCoord) < z_spacings)
                sum += mask(fieldNodes[1][i], fieldNodes[3][j]) * field.data[i, 1, j]
                count += 1*mask(fieldNodes[1][i], fieldNodes[3][j])
            end
        end
    end
    return sum/count
end

#getBackgroundDensity returns the in-situ density of the unaltered water column at depth z using TEOS 10
#The function follows convention such that z is 0 at the surface and becomes increasingly negative with depth. 
function getBackgroundDensity(z, Tfunc::Function, Sfunc::Function)
    eos = TEOS10EquationOfState() 
    S_abs = gsw_sa_from_sp(Sfunc(pipe_radius + 0.0001,z), gsw_p_from_z(z, 30), 31, -30) #31, -30 is a random lat/long in the ocean
    T_consv = gsw_ct_from_t(S_abs, Tfunc(pipe_radius + 0.0001,z), gsw_p_from_z(z, 30)) #set to 30 degrees north
    return TEOS10.ρ(T_consv, S_abs, z, eos) 
end

#checkMemory checks if you have enough RAM to conduct the simulation. it takes in capacity (RAM available), 
#the number of cells in the x-dir, and the number of cells in the z-dir. it will throw an error if you have
#more cells than you have memory for. This is ONLY A ROUGH ESTIMATION, using oceananigans prediction that 
#32gb is approximately 100 million cells. The actual number depends on data saved. 
function checkMemory(capacity, x_cells, z_cells)
    numCells = x_cells * z_cells
    possibleCells = (capacity/32) * 100000000
    if (numCells > possibleCells)
        throw("Too many cells for gpu memory")
        return numCells
    else 
        return numCells
    end
end   



"""CALCULATE DOMAIN DETAILS AND SAVE AS CONSTANTS"""
#OSCILLATION PERIOD DETAILS
#find oscillation period to get relevant timescales & velocities. Approximates via a linear density gradient from top to bottom of pipe
surrounding_density_gradient = (getBackgroundDensity(-pipe_top_depth - pipe_length, TWaterColumn, SWaterColumn) - getBackgroundDensity(-pipe_top_depth, TWaterColumn, SWaterColumn))/pipe_length #takes average for unaltered water column across the pipe (NOT DOMAIN!)
oscillation_angular_frequency =  sqrt((g/1000) * surrounding_density_gradient) #this is N, the buoyancy frequency
oscillation_period = 2π/oscillation_angular_frequency 
@info @sprintf("N = %.3e rad/sec | Buoyancy Oscillation period: %.3f minutes",  oscillation_angular_frequency, oscillation_period/minute)


#VARIOUS GRID SPACING DETAILS 

#option to calculate max grid spacing from predicted maximum velocities 
#max_grid_spacing means its the largest the grid spacing can be, as in, the most coarse the resolution can be. 
#this option below sets it to be just under the decay length scale (4κν/N^2)^(1/4) as derived from balancing momentum and buoyancy equations and assuming exponential velocity decay. Does not seem to be enough to resolve viscous layer. 
#max_grid_spacing = 0.95*((4 * sw_diffusivity_data.T * sw_diffusivity_data.ν)/(oscillation_angular_frequency^2))^(1/4)
@info @sprintf("Max Grid Spacing: %.3em", max_grid_spacing)

#this section sets grid spacing details 
#set to be as big as it can be(to run fastest)
my_x_grid_spacing = max_grid_spacing; #max grid spacing is actually a bit of a misnomer, perhaps shoudl be better called min, its max in the sense that its the max resolution 
my_z_grid_spacing = max_grid_spacing;
#get num of cells in each dimension, ensuring integer number
x_res = floor(Int, domain_x / my_x_grid_spacing);
z_res = floor(Int, domain_z / my_z_grid_spacing);
#width of each grid 
const x_grid_spacing = domain_x/x_res;
const z_grid_spacing = domain_z/z_res;
#call error if it is too many cells to compute 
checkMemory(GPU_memory, x_res, z_res)

#set pipe wall thickness in model, by rounding up to the nearest number of grid cells to desied wall width
const pipe_wall_thickness = roundUp(pipe_wall_thickness_intended, x_grid_spacing)
@info @sprintf("Pipe walls are %.3e meters thick", pipe_wall_thickness)
@info @sprintf("X spacings: %.3e meters | Z spacings: %.3e meters", x_grid_spacing, z_grid_spacing)
@info @sprintf("X resolution: %.3e | Z resolution: %.3e ", x_res, z_res)


#TIME STEPPING DETAILS 
#the goal here is to predict the minimum time step required to fulfill the advective CFL, to allow a relaxation time scale 
#as close to that as possible to reduce leakage of model. Relaxation time scale has to be longer than the longest time step. 

#the option below predicts v max from the buoyancy oscillation, assuming simple harmonic motion. where vmax = angular frequency * amplitude 
#it is no longer used as the parcel is not experience a displaced height as the initial condition. 
#v_max_predicted = oscillation_angular_frequency * height_displaced

#this option below uses a iterative strategy, where you run it once, and then take the max velocity, and use predictions from that
#note that a higher prediced vmax = a shorter max time step (higher temporal resolution). Also note that generally, it seems like you want to 
#use a higher vmax predicted than what is the actual vmax especially at small scales, for higher temporal resolution. 
#0.07 for 20m long, 20cm diameter pipe #0.001 for 1m long 0.02R pipe, 0.001 seems good for small lab scale
v_max_predicted = 0.001 #SIMILITUDE 
min_time_step_predicted = (min(x_grid_spacing, z_grid_spacing)*CFL)/v_max_predicted 
max_time_step_allowed = 2 * min_time_step_predicted #manually set maximum time step  to be 2 * predicted minimum, so that the relaxation time scale can be set close to the time step
#set damping rate for forcers to be 1.1x the maximum time step allowed 
min_relaxation_timescale = 1.1 * max_time_step_allowed
max_relaxation_rate = 1/min_relaxation_timescale


"""MORE COMPUTATIONAL FUNCTIONS THAT ARE DOMAIN DEPENDENT"""
#getXIndex returns the discrete index for x in fieldNodes that is closest to spatial xCoord given
function getXIndex(fieldNodes, xCoord)
    for i in 1 : length(fieldNodes[1])
        if(abs(fieldNodes[1][i] - xCoord) <= 1.01(x_grid_spacing/2)) #1.10 adjusts for some rounding issue with centers 
            return i
        end
    end
    throw("index not found for given z value")
end



"""MASKING FUNCTIONS"""
#=these functions take in spatial coordinates (x, z) or (x, y, z) & return a number between 0 and 1 that describes
how strong the mask is, with 0 being no mask and 1 being full mask. In terms of relaxation, one can think of 
a weaker mask as the equivalent of a longer relaxation timescale, so effectively a more gradual reset to the 
target properties=#

#no mask anywhere
function noMask(x, z)
    return 0
end

#masks out the inside of the center pipe
function isInsidePipe(x, z)
    if (x > (x_center - pipe_radius) && x < (x_center + pipe_radius) && z < -pipe_top_depth && z > -(pipe_top_depth + pipe_length))
        return true
    else
        return false
    end
end
function pipeMask(x, z)
    if (isInsidePipe(x, z))
        return 1
    else
        return 0
    end
end

#masks out the inside of the side pipes 
function isSidePipe(x, z)
    if (min(x, domain_x - x) < pipe_radius && ((-pipe_top_depth) > z > (-pipe_bottom_depth)))
        return true
    else
        return false
    end
end
function sidePipeMask(x, z)
    if (min(x, domain_x - x) < pipe_radius && ((-pipe_top_depth) > z > (-pipe_bottom_depth)))
        return 1
    else
        return 0
    end
end

#option to mask out any interior walls in the middle of the pipe, for example, if you wanted to trigger turbulence 
function isStopper(x, z) #imaginary stopper at center of pipe 
    depth_down = 0.75 # percentage of pipe
    width = 0.6 # percentage of pipe
    thickness = 1 * pipe_wall_thickness
    if ((abs(x_center - x) <= width * pipe_radius) && (abs((-pipe_top_depth - (depth_down * pipe_length)) - z) <= (0.5 * thickness)))
        return true 
    else
        return false 
    end
end

#masks out the pipe walls for center pipe 
function isPipeWall(x, z)
    left_wall_range = (x_center - pipe_radius - pipe_wall_thickness) .. (x_center - pipe_radius)
    right_wall_range = (x_center + pipe_radius) .. (x_center + pipe_radius + pipe_wall_thickness)
    vert_range = -(pipe_top_depth + pipe_length) .. -(pipe_top_depth)
    if ((in(x, left_wall_range) || in(x, right_wall_range)) && in(z, vert_range))
        return true
    #uncomment this following section to incorporate any interior walls 
    # elseif (isStopper(x, z)) 
    #     return true
    else
        return false
    end
end
function pipeWallMask(x,z)
    if (isPipeWall(x, z))
        return 1
    else
        return 0
    end
end

#masks out the pipe walls for the side downwelling pipes. 
function isSidePipeWall(x, z)
    left_wall_range = (pipe_radius) .. (pipe_radius + pipe_wall_thickness)
    right_wall_range = (domain_x - pipe_radius - pipe_wall_thickness) .. (domain_x - pipe_radius)
    vert_range = -(pipe_top_depth + pipe_length) .. -(pipe_top_depth)
    if ((in(x, left_wall_range) || in(x, right_wall_range)) && in(z, vert_range))
        return true
    else
        return false
    end
end
function sidePipeWallMask(x,z)
    if (isSidePipeWall(x, z))
        return 1
    else
        return 0
    end
end

#this mask has the goal of combining all areas that should have a velocity masked to zero. In the current 
#model, this means all areas with a depth between pipe top and pipe bottom that is not inside the center or side pipes. 
function domainVelocityMask(x, z)
    if (!isInsidePipe(x, z) && ((-pipe_bottom_depth) < z < (-pipe_top_depth)) && !isSidePipe(x, z))
        return 1
    elseif (isInsidePipe(x, z) && isPipeWall(x, z)) #this is for the stopper that is for flow sepparation 
        return 1
    else
        return 0
    end
end


#=
These functions describe where in the domain the a field (typically tracer, e.g.) is relaxed to the 
properties of the unaltered water column. 
=#

#tracerRelaxationMaskDomainTwo relaxes all areas below the pipe top that are not a pipe wall or interior
# of a pipe with mask strength 1. Above the pipe, a certain % of the space starting from the pipe
#(as dictated by unmasked ratio is unmasked, and then the the mask goes in a linear gradient from 0 to 0.5
function tracerRelaxationMaskDomainTwo(x, y, z)
    unmasked_ratio = 0.4 # % directly above/below exits of pipe that is not masked 
    if (!isPipeWall(x, y, z) && !isInsidePipe(x, y, z) && (z < (-pipe_top_depth)) && !isSidePipe(x, y, z))
        return 1
    elseif (z > (-pipe_top_depth + unmasked_ratio * pipe_top_depth))
        return (0.5/((1 - unmasked_ratio)*pipe_top_depth))*z + 0.5
    else
        return 0
    end
end

#other method calls
isInsidePipe(x, y, z) = isInsidePipe(x, z)
pipeMask(x,y,z) = pipeMask(x, z)
isSidePipe(x, y, z) = isSidePipe(x, z)
sidePipeMask(x, y, z) = sidePipeMask(x, z)
isPipeWall(x, y, z) = isPipeWall(x, z)
pipeWallMask(x, y, z) = pipeWallMask(x, z)
isSidePipeWall(x, y, z) = isSidePipeWall(x, z)
sidePipeWallMask(x, y, z) = sidePipeWallMask(x, z)
tracerRelaxationMaskDomainTwo(x, z) = tracerRelaxationMaskDomainTwo(x, 0, z)
domainVelocityMask(x, y, z) = domainVelocityMask(x, z)



"""INITIAL CONDITIONS"""
#=these functions take in spatial coordinates (x, z) or (x, y, z) & return the initial conditions of either velocity
or tracer fields at that point.=#


#The goal of neutralDensityPipeTempGradient is to calculate a temperature gradient for the inside of the pipe such that 
#the water in the pipe starts off neutrally buoyant with the surroundings (the density average witin the pipe is the same
#as the density average outside the pipe) However, it does not seem to work. Perhaps since the estimation is based on a linear
#eos? 
#The function takes in the pipe length, the pipe top depth (positive depth), the grid sizing, a function for background water column
#salinity (Sfunc(x,z)) and one for background water column temperature (Tfunc(x, z)) and returns a number for the dT/dz that should
#be in the pipe as the initial condition. 

function neutralDensityPipeTempGradient(p_length, p_top, grid_size, Sfunc::Function, Tfunc::Function)
    eos = TEOS10EquationOfState()
    p_bot = p_top + p_length
    #get reference values at pipe bottom
    ρ₀ = getBackgroundDensity(-p_bot, Tfunc, Sfunc)
    #get average background density & thermal expansion coeff
    z = -p_bot
    density_tot = 0
    count = 0
    α_tot = 0
    while (z <= -p_top)
        density_tot += getBackgroundDensity(z, Tfunc, Sfunc)
        S_abs = gsw_sa_from_sp(Sfunc(0,z), gsw_p_from_z(z, 30), 31, -30) #random lat and long in ocean
        α_tot += SeawaterPolynomials.thermal_expansion(gsw_ct_from_t(S_abs, Tfunc(0, z), gsw_p_from_z(z, 30)), S_abs, z, TEOS10EquationOfState())#set to 30 degrees north
        z += grid_size
        count += 1
    end
         ρ_wc_avg = density_tot/count
    α_avg = α_tot/count   
    #calculate and return the pipe gradient for roughly neutral buoayncy, in the form dt/dz
    return (2 * (ρ₀ - ρ_wc_avg))/(α_avg * 1000 * p_length) #removing the factor of 2 seems to provide a good estimate at least
end 

# const pipeTGrad = neutralDensityPipeTempGradient(pipe_length, pipe_top_depth, max_grid_spacing, SWaterColumn, TWaterColumn)
#unfortunately, due to the failure of my function above, this is currently still set manually :(, equations for testing it are below. 
const pipeTGrad = 0.025#similitude 0.0418 for 1m long

#T_init sets the initial conditions for temperature. Outside of the pipe, this is equal to the water column background. Inside the pipe AND pipe walls
#it follows a gradient (less than the background gradient) such that the pipe (with the low salinity water) is of equivalent densty to the surroundings once 
#integrated through its length
#potentially should change it to just the pipe and not the pipe walls 
function T_init(x, y, z)
    if (isInsidePipe(x, z) || isPipeWall(x, z))
        return TWaterColumn(-pipe_bottom_depth) + (pipeTGrad * (z + pipe_bottom_depth))
    elseif (isSidePipe(x, z))
        return TWaterColumn(-pipe_top_depth) + (pipeTGrad * (z + pipe_top_depth))
    else
        return TWaterColumn(z)
    end
end

#S_init sets the initial conditions for temperature. Outside the pipe and pipe walls, this is equal to the water column background. 
#Inside the pipe AND pipe walls it is equal to the condition at the entrance (low salinity bottom water for the main pipe, high saliity 
#surface water for the side pipes)
function S_init(x, y, z)
    if (isInsidePipe(x, z) || isPipeWall(x, z))
        return SWaterColumn(-pipe_bottom_depth)
    elseif (isSidePipe(x, z))
        return SWaterColumn(-pipe_top_depth)
    else
        return SWaterColumn(z)
    end
end

#A_init sets the initial conditions for arbitrary/visualization tracer A. This sets a value of 1 for the bottom water (all water below pipe 
#bottom, i.e. the water that enters the pipe), and a value of 0 for all else. 
function A_init(x, y, z)
    if(z < -pipe_bottom_depth)
        return 1
    else
        return 0
    end
end 

#w_init sets initial vertical velocity, potentially useful if you wantt to set an initial pumping velocity
function w_init(x, z)
    return 0;
end

#other method calls
T_init(x, z) = T_init(x, 0, z)
S_init(x, z) = S_init(x, 0, z)
A_init(x, z) = A_init(x, 0, z)


#=THE FOLLOWING SECTION IS A DIAGNOSTIC TO CHECK YOUR INITIAL CONDITIONS
Use it to verify the following
    1. the pipe density is similar to the water column density, if it is not exact, try to make the pipe slightly less dense
    this means that wc density - pipe density should be slightly positive, and on the order of 10^-5
    2. wc density - pipe eq density gives you the density of the water column - the density of the fully thermally equillibriated 
    pipe (that is, every point in the pipe has the same temperature as the water at depth z outside the pipe). This is used for 
    similitude testing, so you can use this to check to see if the values are what you expect after salinity gradient scaling, or just
    as another parameter. 
=#

#getDensityPipe returns the density per z coordinate for a x that is located inside the center main pipe (in this case, the very center). 
#this function is used as a diagnostic and takes in a depth z, Tfunc(x, z) and Sfunc(x, z) that return a temperature and salinity respectively
function getDensityPipe(z, Tfunc::Function, Sfunc::Function)
    eos = TEOS10EquationOfState() 
    S_abs = gsw_sa_from_sp(Sfunc(x_center,z), gsw_p_from_z(z, 30), 31, -30) #random lat and long in ocean
    T_consv = gsw_ct_from_t(S_abs, Tfunc(x_center,z), gsw_p_from_z(z, 30)) #set to 30 degrees north
    return TEOS10.ρ(T_consv, S_abs, z, eos) #note that this uses insitu s and T instead of conservative and absolute, which is what the function calls for 
end
z = -pipe_bottom_depth
wc_density_tot = 0
pipe_density_tot = 0
equalized_pipe_density_tot = 0
count = 0
while (z <= -pipe_top_depth)
    wc_density_tot += getBackgroundDensity(z, T_init, S_init) #integrating over non altered water column
    pipe_density_tot += getDensityPipe(z, T_init, S_init) #integrating over pipe initial condition
    equalized_pipe_density_tot += getDensityPipe(z, TWaterColumn, S_init) #integrating over pipe as if fully equalized in temperature
    z += z_grid_spacing
    count += 1
end

@info @sprintf("wc density - pipe_density = %.2e kg/m^3", (wc_density_tot/count) - (pipe_density_tot/count)) # num 1 in the info section above
@info @sprintf("@ eq, wc density - pipe_eq_density = %.3ek kg/m^3", (wc_density_tot/count) - (equalized_pipe_density_tot/count)) #num 2 in the info section above 



"""MODEL BOUNDARY CONDITIONS SETUP"""
#this section calculates the gradients for setting the boundary conditions. As long as you have the base water column functions
#that take in depth input (z), & the z grid spacing, they can be called. 
#TODO: this still need some rigourous testing as I just fixed some syntax confusion
initial_T_top_gradient = (TWaterColumn(0) - TWaterColumn(0 - z_grid_spacing))/z_grid_spacing
initial_T_bottom_gradient = (TWaterColumn(-domain_z + z_grid_spacing) - TWaterColumn(-domain_z))/z_grid_spacing
initial_S_top_gradient = (SWaterColumn(0) - SWaterColumn(0 - z_grid_spacing))/z_grid_spacing
initial_S_bottom_gradient = (SWaterColumn(-domain_z + z_grid_spacing) - SWaterColumn(-domain_z))/z_grid_spacing



"""FORCING FUNCTIONS"""
# This is where you would put forcing functions that you write yourself
# There are a number of example of various forcing functions (discrete and nondiscrete) in the file that was used for 
# prototyping and running on the cpu

#these are diferent method calls, since targets called by relaxation functions pass a (x, z, t) and the water column 
#function takes just a z 
TWaterColumnTarget(x, z, t) = TWaterColumn(z)
SWaterColumnTarget(x, z, t) = SWaterColumn(z)



"""MODEL DIFFUSIVITY SETUP"""
const diffusivity_data = (seawater = sw_diffusivity_data, pipe = pipe_data)
const pipeWallThickness = (actual = pipe_wall_thickness, intended = pipe_wall_thickness_intended)

#= 
The following equations are implemented in the Scalar Diffusivity closure to set different diffusivities and viscosities 
at different spatial locations in the model to allow for changes within the pipe walls. 

In general, they take in the spatial x, y, z coordinates, a named tuple diffusivities which stores the seawater and pipe 
wall material physical data (see above). A number for diffusivity is then returned
=#

#tempDiffusivities returns the thermal diffuivity of the model at each location. Currently, it returned the molecular diffusivity
#at all locations, mimicking an infinite reservoir with water walls, however, it has the capability to return an elevated or reduced
#diffusivity at the pipe walls. 
#tempDiffusivities also takes in wallThickness --> a tuple of information about pipe wall thickness in model vs intended, and 
#a wall indicator which if "WALL" is passed in, the diffusivity for the wall is returned. The wall indicator is used to retrieve
#the wall thermal difffusivity (generally higher) for finding the diffusive CFL later. 
function tempDiffusivities(x, y, z, diffusivities::NamedTuple, wallThickness::NamedTuple, wall_indicator::String)
    return diffusivities[:seawater].T #for infinite reservoir, water walls

    #this following section sets the pipe wall diffusivity to that of the material scaled by the thickness/intended thickness (e.g. if it were
    #twice as thick, it would be 2x the pipe wall material diffusivity)
    # if (isPipeWall(x, z) || wall_indicator == "WALL")
    #     return diffusivities[:pipe].thermal_diffusivity * (wallThickness[:actual]/wallThickness[:intended]) #edit to account for wall thickness
    # else
        #this returns molecular diffusivity
        #return diffusivities[:seawater].T
    # end
end

#saltDiffusivities returns the diffusivity of salt of the model at each location. Currently, this is set to 0 for all pipe walls (side and main)
#and the molecular value in the remainder of the model. 
function saltDiffusivities(x, y, z, diffusivities::NamedTuple)
    if (isPipeWall(x, z) || isSidePipeWall(x, z))
        return 0
    else
        return diffusivities[:seawater].S
    end
end

#myViscosity returns the viscosity of the model at each location. Currently, this is set to 0 for all pipe walls (side and main)
#and the molecular value in the remainder of the model. 
function myViscosity(x, y, z, diffusivities::NamedTuple)
    if (isPipeWall(x, z) || isSidePipeWall(x, z))
        return 0
    else 
        return diffusivities[:seawater].ν
    end

end

#other method calls - these are important, because the ScalarDiffusivity constructor automatically only passes in (x, z, t) where t is time. 
tempDiffusivities(x, z, t) = tempDiffusivities(x, 0, z, diffusivity_data, pipeWallThickness, "")
saltDiffusivities(x ,z ,t ) = saltDiffusivities(x, 0, z, diffusivity_data)
myViscosity(x ,z ,t ) = myViscosity(x, 0, z, diffusivity_data)

"""MODEL SETUP"""
#creates the domain grid. Periodic in x, bounded in z
domain_grid = RectilinearGrid(GPU(), Float64; size=(x_res, z_res), x=(0, domain_x), z=(-domain_z, 0), topology=(Periodic, Flat, Bounded))

#creates clock, starting at time = 0
clock = Clock{eltype(domain_grid)}(time=0)

#timestepper
timestepper = :QuasiAdamsBashforth2

#advection scheme. Default WENO order is 5, can specific higher odd orders. WENO 5 results in slight numerical diffusion 
advection = WENO()

#equation of state. uses TEOS 10
eos = TEOS10EquationOfState()
myBuoyancy = SeawaterBuoyancy(equation_of_state=eos)

#tracers & closure
#MAKE SURE THAT KAPPA IS DEFINED IN THE SAME ORDER AS TRACER DECLARATION
#T is temperature, S is salinity, A is a random tracer for visualization
tracers = (:T, :S, :A)
closure = ScalarDiffusivity(ν=myViscosity, κ=(T=tempDiffusivities, S=saltDiffusivities, A = 0))

#boundary conditions
#sets Neumann (constant gradient) boundary conditions at top/bottom of domain based on initial settings 
T_bcs = FieldBoundaryConditions(top = GradientBoundaryCondition(initial_T_top_gradient), bottom = GradientBoundaryCondition(initial_T_bottom_gradient))
S_bcs = FieldBoundaryConditions(top = GradientBoundaryCondition(initial_S_top_gradient), bottom = GradientBoundaryCondition(initial_S_bottom_gradient))
boundary_conditions = (T = T_bcs, S = S_bcs)

#various forcing functions for right hand side of tracer & momentum equations 
#forces the area covered by domainVelocityMask to have a u and w velocity of 0
uw_domain_forcer = Relaxation(rate = max_relaxation_rate, mask = domainVelocityMask)
#forces the area covered by tracerRelaxationMaskDomainTwo back to background water column temperature properties 
T_domain_forcer = Relaxation(rate = max_relaxation_rate, mask = tracerRelaxationMaskDomainTwo, target = TWaterColumnTarget)
#forces the area covered by tracerRelaxationMaskDomainTwo back to background water column salinity properties 
S_domain_forcer = Relaxation(rate = max_relaxation_rate, mask = tracerRelaxationMaskDomainTwo, target = SWaterColumnTarget)
forcing = (u = (uw_domain_forcer),  w = (uw_domain_forcer), T = (T_domain_forcer), S = (S_domain_forcer))

#instantiates the model
model = NonhydrostaticModel(; grid=domain_grid, clock, advection, buoyancy = myBuoyancy, tracers, timestepper, closure, forcing, boundary_conditions)

#operations to define after model is made. 
density_operation = seawater_density(model; geopotential_height)
@info "model made"

#various checkpoint stuff - implemented but not properly tested!! 
#checkpoint writer 
# model.output_writers[:checkpointer] = Checkpointer(model; frequency = 150000, dir = pathname, prefix = "checkpoint", force = false, verbose = false, properties = [:architecture, :boundary_conditions, :grid, :clock, :buoyancy, :closure, :velocities, :tracers, :timestepper])

#to continue at checkpoint 
# restore_iteration  = 150000
# checkpoint_filename = joinpath(pathname, "checkpoint"*string(restore_iteration)) #TODO: check this 
# model_kwargs = Dict("closure" => closure, "forcing" => forcing)
# model = restore_from_checkpoint(checkpoint_filename; model_kwargs...)



"""SET UP INITIAL CONDITIONS"""
#there was a weird initial bug that required setting the initial conditions twice. this bug is somewhat documented in the github
#effectively, if you define a initial function just before this, you only need to set once. Have not tested in a while - bug may be gone? For now
#seting twice works. 
set!(model, T= T_init, S=S_init, w = w_init, A = A_init)
set!(model, T= T_init, S=S_init, w = w_init, A = A_init)
@info "initial conditions set"



"""SETTING UP SIMULATION"""
#This section calculates various time scales for setting CFL
min_grid_spacing = min(minimum_xspacing(model.grid), minimum_zspacing(model.grid))

#timescales for CFL
viscous_time_scale = (min_grid_spacing^2)/model.closure.ν(0, 0, 0, diffusivity_data)
#for diffusion time scale, use tracer with largest kappa. 
diffusion_time_scale = (min_grid_spacing^2)/model.closure.κ.T(0, 0, 0, diffusivity_data, pipeWallThickness, "WALL") #use for when scalar diffusivity is a function
# diffusion_time_scale = (min_grid_spacing^2)/model.closure.κ.T #use when scalar diffusivity is not a function

#manually set max time step, since diffusion cfl does not work with scalar diffusivity that calls a function to get diffusivity
new_max_time_step = min(0.2 * diffusion_time_scale, 0.2 * viscous_time_scale, max_time_step_allowed) #uses a cfl of 0.2, put in diffusion time scale of tracer with biggest kappa, or viscosity

#timescales for setting initial time step
#the following set it automatically-ish, with just a predicted v_max
# initial_advection_time_scale = min_grid_spacing/v_max_predicted
# initial_time_step = 0.5*min(0.2 * min(diffusion_time_scale, oscillation_period, viscous_time_scale, initial_advection_time_scale), max_time_step_allowed) #sets it to be 1/2 of satisfying all CFLs & max time step allowed
#the following instead allows you to set it manually, which is used now. It shoudl be set to be pretty small to allow the viscous boundary layers to resolve initially. 
initial_time_step = 0.0048 #SIMILITUDE 0.0048 for 1m 0.02R

#initiate simulation
simulation = Simulation(model, Δt=initial_time_step, stop_time=simulation_duration, wall_time_limit=run_duration)



"""SIMULATION CALLBACKS"""
#timewizard for time stepping
timeWizard = TimeStepWizard(cfl=CFL, diffusive_cfl = CFL, max_Δt = new_max_time_step, max_change = 1.01, min_change = 0.5) #max change of 1.01 keeps it more stable
simulation.callbacks[:timeWizard] = Callback(timeWizard, IterationInterval(4))


#progress Message
progress_message(sim) = @printf("Iteration: %04d, time: %s, Δt: %s, wall time: %s\n", iteration(sim), prettytime(sim), prettytime(sim.Δt), prettytime(sim.run_wall_time))
add_callback!(simulation, progress_message, IterationInterval(50))


#output writers

#create model fields 
w = model.velocities.w;
u = model.velocities.u;
T = model.tracers.T;
S = model.tracers.S;
ζ = Field(-∂x(w) + ∂z(u)) #vorticity in y 
ρ = Field(density_operation) #density
A = model.tracers.A

#the following section acts as diagnostics
#it creates fields to store the buoyant and viscous components of the momentum equations for each point in the model.
#note that these are computed manually, not what the model actually calculates, so it may misrepresent the numerics of the model. 
#however, they shoudl serve as a good estimate

#VISCOUS COMPONENT: 
#ViscousComponent returns the portion of a 1D momentum equation that is frictional (viscous) at each discrete point i, j, k. The function takes 
#discrete locations i, j, k, the grid the field is on, the closure of the model, the model clock, and a velocity field (in whichever direction)
#you want the viscosity in. 
@inline function ViscousComponent(i, j, k, grid, closure, clock, velocities)
    #return _νᶜᶜᶜ(i, j, k, grid, model.closure, model.diffusivity_fields, model.clock)  * ∇²ᶜᶜᶜ(i, j, k, grid, model.velocities.w)
    return _νᶜᶜᶜ(i, j, k, grid, closure, nothing, clock)  * ∇²ᶜᶜᶜ(i, j, k, grid, velocities)
end
w_viscousArgs = (model.closure, model.clock, model.velocities.w)
w_viscous_component_kernel_op = KernelFunctionOperation{Center, Center, Center}(ViscousComponent, model.grid, w_viscousArgs...)
ν= Field(w_viscous_component_kernel_op) 
#viscous_field = compute!(ν_component)

#BUOYANCY COMPONENT
#wcBuoyancy returns the portion of the momentum equation (only in w, since the g vector is in just the z) at each discrete point i, j, k. The
#function takes in discrete location i, j, k, the grid the field is on, a SeawaterBuoyancy struct (), and all the tracer fields of the model. Note that this is 
#relative to reference density so it may not be nonzero even if you expect it to be. 
@inline function wcBuoyancy(i, j, k, grid, b::SeawaterBuoyancy, tracer_fields)
    T, S = get_temperature_and_salinity(b, tracer_fields)
    return (- (b.gravitational_acceleration * ρ′(i, j, k, grid, b.equation_of_state, T, S)
    / b.equation_of_state.reference_density))
end

#this computes the buoyancy component in my model relative to the "background density field" - in this case, it subtracts off the density of the water column field
#in an area away from the pipes. As compared to wcIndex, it takes in a 7th argument wc_index, that is the discrete x index of the location to get the "background" from
@inline function BuoyancyComponent(i, j, k, grid, b::SeawaterBuoyancy, tracer_fields, wc_index)
    T, S = get_temperature_and_salinity(b, tracer_fields)
    return (- (b.gravitational_acceleration * ρ′(i, j, k, grid, b.equation_of_state, T, S)
    / b.equation_of_state.reference_density)) - wcBuoyancy(wc_index, j, k, grid, b, tracer_fields)
end
tracer_nodes = nodes(model.grid, (Center(), Center(), Center()), reshape = true)
wc_index = getXIndex(tracer_nodes, pipe_radius + (3 *x_grid_spacing))

#this is the nonlinear option that subtracts off the background
buoyancyArgs = (myBuoyancy, model.tracers, wc_index) #for some reason model.buoyancy.equation_of_state does not return an eos, so I use myBuoyancy

#this is the built in buoyancy perturbation option that only takes a linear eos (otherwise identicle to wcBuoyancy)
# buoyancyArgs_linear = (SeawaterBuoyancy(), model.tracers)
# #buoyancy_component_kernel_op = KernelFunctionOperation{Center, Center, Center}(buoyancy_perturbationᶜᶜᶜ, model.grid, buoyancyArgs_linear...) #uses built in function which is linear

buoyancy_component_kernel_op = KernelFunctionOperation{Center, Center, Center}(BuoyancyComponent, model.grid, buoyancyArgs...)
b = Field(buoyancy_component_kernel_op)
#buoyancy_field = compute!(b_component)

#declare output writers
#option with diagnostics
simulation.output_writers[:outputs] = JLD2OutputWriter(model, (; u, w, T, S, ζ, ρ, A, ν, b); filename, schedule=TimeInterval(output_interval), overwrite_existing=true) #can also set to TimeInterval
#option without diagnostics 
#simulation.output_writers[:outputs] = JLD2OutputWriter(model, (; u, w, T, S, ζ, ρ, A); filename, schedule=TimeInterval(output_interval), overwrite_existing=true) #can also set to TimeInterval

#run simulation
run!(simulation, pickup = false)

#ability to keep running, note that it starts from wall time 0 again, but simulation time is whatever it left off as 
# simulation.wall_time_limit = 10hour
# run!(simulation)




"""DATA RETREIVAL FROM FILE"""
output_filename = filename * ".jld2"
u_t = FieldTimeSeries(output_filename, "u")
w_t = FieldTimeSeries(output_filename, "w")
ζ_t = FieldTimeSeries(output_filename, "ζ")
T_t = FieldTimeSeries(output_filename, "T")
S_t = FieldTimeSeries(output_filename, "S")
ρ_t = FieldTimeSeries(output_filename, "ρ")
A_t = FieldTimeSeries(output_filename, "A")
b_t = FieldTimeSeries(output_filename, "b")
ν_t = FieldTimeSeries(output_filename, "ν")
times = u_t.times
n = Observable(1)
Observable(1)
@info @sprintf("Finished extracting time series. Simulation total run time %.3f seconds | %.3f minutes", times[end], times[end]/minute)
uₙ = @lift interior(u_t[$n], :, 1, :)
wₙ = @lift interior(w_t[$n], :, 1, :)
ζₙ = @lift interior(ζ_t[$n], :, 1, :)
Tₙ = @lift interior(T_t[$n], :, 1, :)
Sₙ = @lift interior(S_t[$n], :, 1, :)
ρₙ = @lift interior(ρ_t[$n], :, 1, :)
Aₙ = @lift interior(A_t[$n], :, 1, :)
bₙ = @lift interior(b_t[$n], :, 1, :)
νₙ = @lift interior(ν_t[$n], :, 1, :)
bνₙ= @lift (interior(ν_t[$n], :, 1, :) + interior(b_t[$n], :, 1, :))
#how much of data set to plot 
num_Data_Points = length(times)

#very inefficient way of getting max/min, need to update
T_range = getMaxAndMinEnd(T_t)
S_range = getMaxAndMinEnd(S_t)
ρ_range = getMaxAndMinEnd(ρ_t)
u_range = getMaxAndMinEnd(u_t)
w_range = getMaxAndMinEnd(w_t)
ζ_range = getMaxAndMinEnd(ζ_t)
A_range = getMaxAndMinEnd(A_t)
# T_range = getMaxAndMin(num_Data_Points, T_t)
# S_range = getMaxAndMin(num_Data_Points, S_t)
# ρ_range = getMaxAndMin(num_Data_Points, ρ_t)
# u_range = getMaxAndMin(num_Data_Points, u_t)
# w_range = getMaxAndMin(num_Data_Points, w_t)
# ζ_range = getMaxAndMin(num_Data_Points, ζ_t)
# A_range = getMaxAndMin(num_Data_Points, A_t)
@info "finished getting max and min of each"



""""PLOTTING DATA"""
#plotting the average along the pipe vs time 
function plotAverages(timeSeriesField, myTimes, locs)
    averages = zeros(length(myTimes))
    colors = (:lightblue, :blue, :navy)
    for i in 0.25:0.25:0.75
        for j in 1 : length(myTimes)
            averages[j] = getMaskedAverageAtZ(pipeMask, timeSeriesField[j], locs, (-pipe_top_depth - i*pipe_length), z_grid_spacing)
        end
        scatterlines!((myTimes/hour), averages, label = @sprintf("@ %1.2f pipe", i), color = colors[i/0.25], markersize = 5)
    end

    average_all = zeros(length(myTimes))
    for i in 1 : length(myTimes)
        average_all[i] = getMaskedAverage(pipeMask, timeSeriesField[i], locs)
    end
    scatterlines!((myTimes/hour), average_all, label = "avg entire pipe", color = :red, marker = :cross, markersize = 5)
end

#AVERAGE VALUES
fig = Figure(size = (1000, 600))
#fig[0, :] = Label(fig, title)
title = "Averaged Values"
fig[0, :] = Label(fig, title)
velocities_plot= Axis(fig[1,1], title = "Pipe Velocities Averaged", xlabel="time(hrs)", ylabel = "velocities (m/s)", width =  700)
plotAverages(w_t, times, (Center(), Center(), Face()))
xlims!(0, times[end]/hour)
fig[1, 2] = Legend(fig, velocities_plot, frame_visible = false)
temps_plot = Axis(fig[2,1], title = "Pipe Temperatures Averaged", xlabel="time(hrs)", ylabel = "Temperature(C)", width =  700)
plotAverages(T_t, times, (Center(), Center(), Center()))
xlims!(0, times[end]/hour)
scatterlines!(times, fill(TWaterColumn(0), length(times)), label = "surface value", color = :maroon1, markersize = 0, linestyle = :dashdotdot,)
scatterlines!(times, fill(TWaterColumn(-domain_z), length(times)), label = "bottom value", color = :mediumpurple2, markersize = 0, linestyle = :dashdotdot)
fig[2, 2] = Legend(fig, temps_plot, frame_visible = false)
salinities_plot = Axis(fig[3, 1], title = "Pipe Salinity Averaged", xlabel="time(hrs)", ylabel = "Salinity(ppt)", width =  700)
plotAverages(S_t, times, (Center(), Center(), Center()))
xlims!(0, times[end]/hour)
scatterlines!(times, fill(SWaterColumn(0), length(times)), label = "surface value", color = :maroon1, markersize = 0, linestyle = :dashdotdot)
scatterlines!(times, fill(SWaterColumn(-domain_z), length(times)), label = "bottom value", color = :mediumpurple2, markersize = 0, linestyle = :dashdotdot)
fig[3, 2] = Legend(fig, salinities_plot, frame_visible = false)
fig
save(joinpath(pathname,"averaged values vs time.png"), fig)
@info "Finished plotting average values vs time chart"


"""PLOTTING CROSS SECTIONAL STUFF VS TIME""" 
#returns indexes [first index of left wall, first index of pipe, last index of pipe, last index of right wall]
function getPipeAndWallXIndexRange(fieldNodes)
    leftWallIndex = -1
    rightWallIndex = -1
    leftPipeIndex = -1
    rightPipeIndex = -1
    z_center_index = ceil(Int, length(fieldNodes[3])/2)
    edgedetector = 0;
    count = 0
    for i in 1 : length(fieldNodes[1])
        if (pipeWallMask(fieldNodes[1][i], fieldNodes[3][z_center_index]) != edgedetector)
            if (edgedetector == 0)
                count += 1
                if (count == 1)
                    leftWallIndex = i
                else
                    rightPipeIndex = i - 1
                end
                edgedetector = 1
            else 
                count += 1
                if (count == 4)
                    rightWallIndex = i - 1
                    break;
                else
                    leftPipeIndex = i 
                end
                edgedetector = 0
            end
        end
    end
    if (leftWallIndex == -1 || rightWallIndex == -1 || leftPipeIndex == -1 || rightPipeIndex == -1)
        throw("pipe side walls not found")
    else
        return [leftWallIndex, leftPipeIndex, rightPipeIndex, rightWallIndex]
    end
end 
#This functon returns the discrete index for z, rounding to the nearest one. If the point is perfectly in between two, it returns the lower one 
function getZIndex(fieldNodes, zCoord)
    for i in 1 : length(fieldNodes[3])
        if(abs(fieldNodes[3][i] - zCoord) <= 1.01(z_grid_spacing/2)) #1.10 adjusts for some rounding issue with centers 
            return i
        end
    end
    throw("index not found for given z value")
end

#general info
w_velocity_nodes = nodes(w_t[1].grid, (Center(), Center(), Face()), reshape = true)
x_pipe_range_velocities= getPipeAndWallXIndexRange(w_velocity_nodes) 
x_plot_range_velocities = w_velocity_nodes[1][(x_pipe_range_velocities[1] - 4):(x_pipe_range_velocities[4] + 4)]#gives a 4 cell halo
z_index_tenth_velocities = getZIndex(w_velocity_nodes, -pipe_top_depth - (0.1 * pipe_length)) 
z_index_quarter_velocities = getZIndex(w_velocity_nodes, -pipe_top_depth - (0.25 * pipe_length)) 
z_index_half_velocities = getZIndex(w_velocity_nodes, -pipe_top_depth - (0.5 * pipe_length)) 
z_index_three_quarter_velocities = getZIndex(w_velocity_nodes, -pipe_top_depth - (0.75 * pipe_length)) 
tracer_nodes = nodes(T_t[1].grid, (Center(), Center(), Center()), reshape = true)
x_pipe_range_tracer= getPipeAndWallXIndexRange(tracer_nodes) 
x_plot_range_tracer = tracer_nodes[1][(x_pipe_range_tracer[1] - 4):(x_pipe_range_tracer[4] + 4)]#gives a 4 cell halo
z_index_tenth_tracer = getZIndex(tracer_nodes, -pipe_top_depth - (0.1 * pipe_length)) 
z_index_quarter_tracer = getZIndex(tracer_nodes, -pipe_top_depth - (0.25 * pipe_length)) 
z_index_half_tracer = getZIndex(tracer_nodes, -pipe_top_depth - (0.5 * pipe_length)) 
z_index_three_quarter_tracer = getZIndex(tracer_nodes, -pipe_top_depth - (0.75 * pipe_length)) 

#extract cross sections
wₙ_tenth =  @lift interior(w_t[$n], (x_pipe_range_velocities[1] - 4):(x_pipe_range_velocities[4] + 4) , 1, z_index_tenth_velocities)
wₙ_quarter =  @lift interior(w_t[$n], (x_pipe_range_velocities[1] - 4):(x_pipe_range_velocities[4] + 4) , 1, z_index_quarter_velocities)
wₙ_half =  @lift interior(w_t[$n], (x_pipe_range_velocities[1] - 4):(x_pipe_range_velocities[4] + 4) , 1, z_index_half_velocities)
wₙ_three_quarter =  @lift interior(w_t[$n], (x_pipe_range_velocities[1] - 4):(x_pipe_range_velocities[4] + 4) , 1, z_index_three_quarter_velocities)

#get max and min in pipe
function getMaxAndMinInPipe(numPoints, dataSeries, x_range::UnitRange, z_range::UnitRange)
    myMax = maximum(interior(dataSeries[1], x_range, 1, z_range))
    myMin = minimum(interior(dataSeries[1], x_range, 1, z_range))
    for i in 2:numPoints
        myMax = max(maximum(interior(dataSeries[i], x_range, 1, z_range)),myMax)
        myMin = min(minimum(interior(dataSeries[i], x_range, 1, z_range)),myMin)
    end
    return (myMin, myMax)
end
w_pipe_range = getMaxAndMinInPipe(num_Data_Points, w_t, (x_pipe_range_velocities[2]):(x_pipe_range_velocities[3]), (getZIndex(w_velocity_nodes, -pipe_bottom_depth + 0.01)):(getZIndex(w_velocity_nodes, -pipe_top_depth - 0.01)))

fig = Figure(size = (1000, 600))
# myTicksBIG = -pipe_radius : pipe_radius/10: pipe_radius
# myTicksSMALL = -pipe_radius : pipe_radius/5: pipe_radius
kwargsBIG = (; xminorticks = IntervalsBetween(10), xminorticksvisible = true, xminorgridvisible = true)
kwargs =(; xminorticks = IntervalsBetween(5), xminorticksvisible = true, xminorgridvisible = true)
title = @lift @sprintf("t = %1.2f minutes", round(times[$n] / minute, digits=2))
fig[0, 1:4] = Label(fig, title)
cross_velocities_plot= Axis(fig[1,1:3], title = "Pipe Cross Sectional Velocities", xlabel="x (m), centered at pipe center", ylabel = "velocities (m/s)", width =  600; kwargsBIG...)
scatterlines!((x_plot_range_velocities .- x_center), wₙ_tenth, label = "@ 0.1 pipe", color = colors = :gray, markersize = 5)
scatterlines!((x_plot_range_velocities .- x_center), wₙ_quarter, label = "@ 0.25 pipe", color = colors = :lightblue, markersize = 5)
scatterlines!((x_plot_range_velocities .- x_center), wₙ_half, label = "@ 0.5 pipe", color = colors = :blue, markersize = 5)
scatterlines!((x_plot_range_velocities .- x_center), wₙ_three_quarter, label = "@ 0.75 pipe", color = colors = :navy, markersize = 5)
vlines!([(w_velocity_nodes[1][x_pipe_range_velocities[4]] - x_center),(w_velocity_nodes[1][x_pipe_range_velocities[1]] - x_center),(w_velocity_nodes[1][x_pipe_range_velocities[2] - 1] - x_center) , (w_velocity_nodes[1][x_pipe_range_velocities[3] + 1] - x_center)], label = "pipe side walls", color = :deeppink2, linewidth = 5)
ylims!(-max(abs(w_pipe_range[1]), abs(w_pipe_range[2])), max(abs(w_pipe_range[1]), abs(w_pipe_range[2])))
# ylims!(-0.01, 0.01)
fig[1, 4] = Legend(fig, cross_velocities_plot, frame_visible = false)
quarterPlot= Axis(fig[2, 1], width = 150, title = "@ 0.1 pipe"; kwargs...)
scatterlines!((x_plot_range_velocities .- x_center), wₙ_tenth, color = :gray, markersize = 3)
vlines!([(w_velocity_nodes[1][x_pipe_range_velocities[4]] - x_center),(w_velocity_nodes[1][x_pipe_range_velocities[1]] - x_center),(w_velocity_nodes[1][x_pipe_range_velocities[2] - 1] - x_center) , (w_velocity_nodes[1][x_pipe_range_velocities[3] + 1] - x_center)], label = "pipe side walls", color = :deeppink2, linewidth = 3)
ylims!(-max(abs(w_pipe_range[1]), abs(w_pipe_range[2])), max(abs(w_pipe_range[1]), abs(w_pipe_range[2])))
# ylims!(-0.01, 0.01)
quarterPlot= Axis(fig[2, 2], width = 150, title = "@ 0.25 pipe"; kwargs...)
scatterlines!((x_plot_range_velocities .- x_center), wₙ_quarter, color = :lightblue, markersize = 3)
vlines!([(w_velocity_nodes[1][x_pipe_range_velocities[4]] - x_center),(w_velocity_nodes[1][x_pipe_range_velocities[1]] - x_center),(w_velocity_nodes[1][x_pipe_range_velocities[2] - 1] - x_center) , (w_velocity_nodes[1][x_pipe_range_velocities[3] + 1] - x_center)], label = "pipe side walls", color = :deeppink2, linewidth = 3)
ylims!(-max(abs(w_pipe_range[1]), abs(w_pipe_range[2])), max(abs(w_pipe_range[1]), abs(w_pipe_range[2])))
# ylims!(-0.01, 0.01)
quarterPlot= Axis(fig[2, 3], width = 150, title = "@ 0.5 pipe"; kwargs...)
scatterlines!((x_plot_range_velocities .- x_center), wₙ_half, color = :blue, markersize = 3)
vlines!([(w_velocity_nodes[1][x_pipe_range_velocities[4]] - x_center),(w_velocity_nodes[1][x_pipe_range_velocities[1]] - x_center),(w_velocity_nodes[1][x_pipe_range_velocities[2] - 1] - x_center) , (w_velocity_nodes[1][x_pipe_range_velocities[3] + 1] - x_center)], label = "pipe side walls", color = :deeppink2, linewidth = 3)
ylims!(-max(abs(w_pipe_range[1]), abs(w_pipe_range[2])), max(abs(w_pipe_range[1]), abs(w_pipe_range[2])))
# ylims!(-0.01, 0.01)
quarterPlot= Axis(fig[2, 4], width = 150, title = "@ 0.75 pipe"; kwargs...)
scatterlines!((x_plot_range_velocities .- x_center), wₙ_three_quarter, color = :navy, markersize = 3)
vlines!([(w_velocity_nodes[1][x_pipe_range_velocities[4]] - x_center),(w_velocity_nodes[1][x_pipe_range_velocities[1]] - x_center),(w_velocity_nodes[1][x_pipe_range_velocities[2] - 1] - x_center) , (w_velocity_nodes[1][x_pipe_range_velocities[3] + 1] - x_center)], label = "pipe side walls", color = :deeppink2, linewidth = 3)
ylims!(-max(abs(w_pipe_range[1]), abs(w_pipe_range[2])), max(abs(w_pipe_range[1]), abs(w_pipe_range[2])))
# ylims!(-0.01, 0.01)
fig
@info "Making cross sectional velocties animation from data"
frames = 1:length(times)
record(fig, joinpath(pathname,"CROSSSECvelocities.mp4"), frames, framerate=8) do i
    @info string("Plotting frame ", i, " of ", frames[end])
    n[] = i
end


#temperatures
#extract cross sections
Tₙ_quarter =  @lift interior(T_t[$n], (x_pipe_range_tracer[1] - 4):(x_pipe_range_tracer[4] + 4) , 1, z_index_quarter_tracer)
Tₙ_half =  @lift interior(T_t[$n], (x_pipe_range_tracer[1] - 4):(x_pipe_range_tracer[4] + 4) , 1, z_index_half_tracer)
Tₙ_three_quarter =  @lift interior(T_t[$n], (x_pipe_range_tracer[1] - 4):(x_pipe_range_tracer[4] + 4) , 1, z_index_three_quarter_tracer)

fig = Figure(size = (700, 600))
title = @lift @sprintf("t = %1.2f minutes", round(times[$n] / minute, digits=2))
fig[0, 1:3] = Label(fig, title)
cross_tracer_plot= Axis(fig[1,1:2], title = "Pipe Cross Sectional temperatures", xlabel="x (m), centered at pipe center", ylabel = "Temperature (C)", width =  400)
scatterlines!((x_plot_range_tracer .- x_center), Tₙ_quarter, label = "@ 0.25 pipe", color = colors = :lightblue, markersize = 5)
scatterlines!((x_plot_range_tracer .- x_center), Tₙ_half, label = "@ 0.5 pipe", color = colors = :blue, markersize = 5)
scatterlines!((x_plot_range_tracer .- x_center), Tₙ_three_quarter, label = "@ 0.75 pipe", color = colors = :navy, markersize = 5)
vlines!([(w_velocity_nodes[1][x_pipe_range_tracer[4]] - x_center),(w_velocity_nodes[1][x_pipe_range_tracer[1]] - x_center),(w_velocity_nodes[1][x_pipe_range_tracer[2] - 1] - x_center) , (w_velocity_nodes[1][x_pipe_range_tracer[3] + 1] - x_center)], label = "pipe side walls", color = :deeppink2, linewidth = 5)
ylims!(T_range[1], T_range[2])
fig[1, 3] = Legend(fig, cross_tracer_plot, frame_visible = false)
quarterPlot= Axis(fig[2, 1], width = 150, title = "@ 0.25 pipe")
scatterlines!((x_plot_range_tracer .- x_center), Tₙ_quarter, color = :lightblue, markersize = 3)
vlines!([(w_velocity_nodes[1][x_pipe_range_tracer[4]] - x_center),(w_velocity_nodes[1][x_pipe_range_tracer[1]] - x_center),(w_velocity_nodes[1][x_pipe_range_tracer[2] - 1] - x_center) , (w_velocity_nodes[1][x_pipe_range_tracer[3] + 1] - x_center)], label = "pipe side walls", color = :deeppink2, linewidth = 3)
ylims!(T_range[1], T_range[2])
quarterPlot= Axis(fig[2, 2], width = 150, title = "@ 0.5 pipe")
scatterlines!((x_plot_range_tracer .- x_center), Tₙ_half, color = :blue, markersize = 3)
vlines!([(w_velocity_nodes[1][x_pipe_range_tracer[4]] - x_center),(w_velocity_nodes[1][x_pipe_range_tracer[1]] - x_center),(w_velocity_nodes[1][x_pipe_range_tracer[2] - 1] - x_center) , (w_velocity_nodes[1][x_pipe_range_tracer[3] + 1] - x_center)], label = "pipe side walls", color = :deeppink2, linewidth = 3)
ylims!(T_range[1], T_range[2])
quarterPlot= Axis(fig[2, 3], width = 150, title = "@ 0.75 pipe")
scatterlines!((x_plot_range_tracer .- x_center), Tₙ_three_quarter, color = :navy, markersize = 3)
vlines!([(w_velocity_nodes[1][x_pipe_range_tracer[4]] - x_center),(w_velocity_nodes[1][x_pipe_range_tracer[1]] - x_center),(w_velocity_nodes[1][x_pipe_range_tracer[2] - 1] - x_center) , (w_velocity_nodes[1][x_pipe_range_tracer[3] + 1] - x_center)], label = "pipe side walls", color = :deeppink2, linewidth = 3)
ylims!(T_range[1], T_range[2])
fig
@info "Making cross sectional temperature animation from data"
frames = 1:length(times)
record(fig, joinpath(pathname, "EXTENDEDCROSSSECtemperatures.mp4"), frames, framerate=8) do i
    @info string("Plotting frame ", i, " of ", frames[end])
    n[] = i
end



"""STANDARD ANIMATIONS"""
#properties
fig = Figure(size=(600, 900))
title = @lift @sprintf("t = %1.2f minutes", round(times[$n] / minute, digits=2))
axis_kwargs = (xlabel="x (m)", ylabel="z (m)", width=400)
fig[1, :] = Label(fig, title)
xT, yT, zT = nodes(T_t[1])
ax_T = Axis(fig[2, 1]; title="temperature[C]", axis_kwargs...)
hm_T = heatmap!(ax_T, xT, zT, Tₙ; colorrange=T_range, colormap=:thermal) #note that this is still using old grid from T, S, initial, may need to recompute x and z using specific nodes 
Colorbar(fig[2, 2], hm_T, label="C")
xS, yS, zS = nodes(S_t[1])
ax_S = Axis(fig[3, 1]; title="salinity[ppt]", axis_kwargs...)
hm_S = heatmap!(ax_S, xS, zS, Sₙ; colorrange=S_range, colormap=:haline) #note that this is still using old grid from T, S, initial, may need to recompute x and z using specific nodes 
Colorbar(fig[3, 2], hm_S, label="ppt")
xρ, yρ, zρ = nodes(ρ_t[1])
ax_ρ = Axis(fig[4, 1]; title="potential density[kg/m^3]", axis_kwargs...)
hm_ρ = heatmap!(ax_ρ, xρ, zρ, ρₙ; colorrange=ρ_range, colormap=Reverse(:viridis)) #note that this is still using old grid from T, S, initial, may need to recompute x and z using specific nodes 
Colorbar(fig[4, 2], hm_ρ, label="kg/m^3")
fig
@info "Making properties animation from data"
frames = 1:length(times)
record(fig, joinpath(pathname,"properties.mp4"), frames, framerate=8) do i
    @info string("Plotting frame ", i, " of ", frames[end])
    n[] = i
end

#velocities
fig = Figure(size=(600, 900))
title = @lift @sprintf("t = %1.2f minutes", round(times[$n] / minute, digits=2))
axis_kwargs = (xlabel="x (m)", ylabel="z (m)", width=400)
fig[1, :] = Label(fig, title)
xw, yw, zw = nodes(w_t[1])
ax_w = Axis(fig[2, 1]; title="w velocity", axis_kwargs...)
w_colorbar_range = (-max(abs(w_range[1]), abs(w_range[2])), max(abs(w_range[1]), abs(w_range[2])))
hm_w = heatmap!(ax_w, xw, zw, wₙ; colorrange=w_colorbar_range, colormap=:balance) #note that this is still using old grid from T, S, initial, may need to recompute x and z using specific nodes 
Colorbar(fig[2, 2], hm_w, label="m/s")
xu, yu, zu = nodes(u_t[1])
ax_u = Axis(fig[3, 1]; title="u velocity", axis_kwargs...)
u_colorbar_range = (-max(abs(u_range[1]), abs(u_range[2])), max(abs(u_range[1]), abs(u_range[2])))
hm_u = heatmap!(ax_u, xu, zu, uₙ; colorrange=u_colorbar_range, colormap=:balance) #note that this is still using old grid from T, S, initial, may need to recompute x and z using specific nodes 
Colorbar(fig[3, 2], hm_u, label="m/s")
xζ, yζ, zζ = nodes(ζ_t[1])
ax_ζ = Axis(fig[4, 1]; title="vorticity", axis_kwargs...)
ζ_colorbar_range = (-max(abs(ζ_range[1]), abs(ζ_range[2])), max(abs(ζ_range[1]), abs(ζ_range[2])))
hm_ζ = heatmap!(ax_ζ, xζ, zζ, ζₙ; colorrange=ζ_colorbar_range, colormap=:balance) #note that this is still using old grid from T, S, initial, may need to recompute x and z using specific nodes 
Colorbar(fig[4, 2], hm_ζ, label="rot/sec")
fig
@info "Making velocities animation from data"
frames = 1:length(times)
record(fig, joinpath(pathname,"velocities.mp4"), frames, framerate=8) do i
    @info string("Plotting frame ", i, " of ", frames[end])
    n[] = i
end

#other tracer
fig = Figure(size=(600, 300))
title = @lift @sprintf("t = %1.2f minutes", round(times[$n] / minute, digits=2))
axis_kwargs = (xlabel="x (m)", ylabel="z (m)", width=400)
fig[1, :] = Label(fig, title)
xA, yA, zA = nodes(A_t[1])
ax_A = Axis(fig[2, 1]; title="tracer A", axis_kwargs...)
A_colorbar_range = (A_range)
hm_A = heatmap!(ax_A, xA, zA, Aₙ; colorrange=A_colorbar_range, colormap=:matter) 
Colorbar(fig[2, 2], hm_A, label="amount")
frames = 1:length(times)
@info "Making tracers animation from data"
record(fig, joinpath(pathname,"tracer.mp4"), frames, framerate=8) do i
    @info string("Plotting frame ", i, " of ", frames[end])
    n[] = i
end
@info @sprintf("Done. Simulation total run time %.3f seconds | %.3f minutes", times[end], times[end]/minute)




"""PLOTTING ADVANCED (ZOOMED)"""
# just ploting area around pipe with 5x pipe diameteron either side & 2x pipe diameter on top and bottom
#setting ranges, for testing, should delete 
#velocities
fig = Figure(size=(900, 600))
title = @lift @sprintf("t = %1.2f minutes", round(times[$n] / minute, digits=2))
axis_kwargs = (xlabel="x (m)", ylabel="z (m)", width=300)
fig[1, 1:3] = Label(fig, title)
xw, yw, zw = nodes(w_t[1])
ax_w = Axis(fig[2, 1]; title="w velocity", axis_kwargs...)
w_colorbar_range = (-max(abs(w_range[1]), abs(w_range[2])), max(abs(w_range[1]), abs(w_range[2])))
hm_w = heatmap!(ax_w, xw, zw, wₙ; colorrange=w_colorbar_range, colormap=:balance) #note that this is still using old grid from T, S, initial, may need to recompute x and z using specific nodes 
xlims!(ax_w,(x_center - pipe_radius - pipe_wall_thickness - (1.5 * pipe_radius)), (x_center + pipe_radius + pipe_wall_thickness + (1.5 * pipe_radius)))
ylims!(ax_w, (-pipe_bottom_depth - (8 * pipe_radius)), (-pipe_top_depth + (8 * pipe_radius))) 
vlines!([(w_velocity_nodes[1][x_pipe_range_velocities[4]]),(w_velocity_nodes[1][x_pipe_range_velocities[1]]),(w_velocity_nodes[1][x_pipe_range_velocities[2] - 1]) , (w_velocity_nodes[1][x_pipe_range_velocities[3] + 1])], label = "pipe side walls", color = :deeppink2, linewidth = 1)
#Colorbar(fig[2, 2], hm_w, label="m/s")
xu, yu, zu = nodes(u_t[1])
ax_u = Axis(fig[2, 2]; title="u velocity", axis_kwargs...)
u_colorbar_range = (-max(abs(u_range[1]), abs(u_range[2])), max(abs(u_range[1]), abs(u_range[2])))
hm_u = heatmap!(ax_u, xu, zu, uₙ; colorrange=w_colorbar_range, colormap=:balance) #note that this is still using old grid from T, S, initial, may need to recompute x and z using specific nodes 
xlims!(ax_u,(x_center - pipe_radius - pipe_wall_thickness - (1.5 * pipe_radius)), (x_center + pipe_radius + pipe_wall_thickness + (1.5 * pipe_radius)))
ylims!(ax_u, (-pipe_bottom_depth - (8 * pipe_radius)), (-pipe_top_depth + (8 * pipe_radius))) 
vlines!([(w_velocity_nodes[1][x_pipe_range_velocities[4]]),(w_velocity_nodes[1][x_pipe_range_velocities[1]]),(w_velocity_nodes[1][x_pipe_range_velocities[2] - 1]) , (w_velocity_nodes[1][x_pipe_range_velocities[3] + 1])], label = "pipe side walls", color = :deeppink2, linewidth = 1)
Colorbar(fig[2, 3], hm_u, label="m/s")
# xζ, yζ, zζ = nodes(ζ_t[1])
# ax_ζ = Axis(fig[4, 1]; title="vorticity", axis_kwargs...)
# ζ_colorbar_range = (-max(abs(ζ_range[1]), abs(ζ_range[2])), max(abs(ζ_range[1]), abs(ζ_range[2])))
# hm_ζ = heatmap!(ax_ζ, xζ, zζ, ζₙ; colorrange=ζ_colorbar_range, colormap=:balance) #note that this is still using old grid from T, S, initial, may need to recompute x and z using specific nodes 
# xlims!(ax_ζ,(x_center - pipe_radius - pipe_wall_thickness - (5 * pipe_radius)), (x_center + pipe_radius + pipe_wall_thickness + (5 * pipe_radius)))
# ylims!(ax_ζ, (-pipe_bottom_depth - (8 * pipe_radius)), (-pipe_top_depth + (8 * pipe_radius))) 
# vlines!([(w_velocity_nodes[1][x_pipe_range_velocities[4]]),(w_velocity_nodes[1][x_pipe_range_velocities[1]]),(w_velocity_nodes[1][x_pipe_range_velocities[2] - 1]) , (w_velocity_nodes[1][x_pipe_range_velocities[3] + 1])], label = "pipe side walls", color = :deeppink2, linewidth = 1)
# Colorbar(fig[4, 2], hm_ζ, label="rot/sec")
fig
@info "Making velocities animation from data"
frames = 1:length(times)
record(fig, joinpath(pathname,"ZOOMvelocities.mp4"), frames, framerate=8) do i
    @info string("Plotting frame ", i, " of ", frames[end])
    n[] = i
end

#tracers
fig = Figure(size=(600, 300))
title = @lift @sprintf("t = %1.2f minutes", round(times[$n] / minute, digits=2))
axis_kwargs = (xlabel="x (m)", ylabel="z (m)", width=150)
fig[1, :] = Label(fig, title)
xA, yA, zA = nodes(A_t[1])
ax_A = Axis(fig[2, 1]; title="tracer A", axis_kwargs...)
A_colorbar_range = (A_range)
hm_A = heatmap!(ax_A, xA, zA, Aₙ; colorrange=A_colorbar_range, colormap=:matter) 
xlims!(ax_A,(x_center - pipe_radius - pipe_wall_thickness - (5 * pipe_radius)), (x_center + pipe_radius + pipe_wall_thickness + (5 * pipe_radius)))
ylims!(ax_A, (-pipe_bottom_depth - (8 * pipe_radius)), (-pipe_top_depth + (8 * pipe_radius))) 
vlines!([(w_velocity_nodes[1][x_pipe_range_velocities[4]]),(w_velocity_nodes[1][x_pipe_range_velocities[1]]),(w_velocity_nodes[1][x_pipe_range_velocities[2] - 1]) , (w_velocity_nodes[1][x_pipe_range_velocities[3] + 1])], label = "pipe side walls", color = :deeppink2, linewidth = 1)
Colorbar(fig[2, 2], hm_A, label="amount")
frames = 1:length(times)
@info "Making tracers animation from data"
record(fig, joinpath(pathname,"ZOOMtracer.mp4"), frames, framerate=8) do i
    @info string("Plotting frame ", i, " of ", frames[end])
    n[] = i
end

#properties
fig = Figure(size=(1300, 600))
title = @lift @sprintf("t = %1.2f minutes", round(times[$n] / minute, digits=2))
axis_kwargs = (xlabel="x (m)", width=300)
fig[1, 1:6] = Label(fig, title)
xT, yT, zT = nodes(T_t[1])
ax_T = Axis(fig[2, 1]; title="temperature[C]", ylabel="z (m)", axis_kwargs...)
hm_T = heatmap!(ax_T, xT, zT, Tₙ; colorrange=T_range, colormap=:thermal)
xlims!(ax_T,(x_center - pipe_radius - pipe_wall_thickness - (1.5 * pipe_radius)), (x_center + pipe_radius + pipe_wall_thickness + (1.5 * pipe_radius)))
ylims!(ax_T, (-pipe_bottom_depth - (8 * pipe_radius)), (-pipe_top_depth + (8 * pipe_radius)))  #note that this is still using old grid from T, S, initial, may need to recompute x and z using specific nodes 
vlines!([(w_velocity_nodes[1][x_pipe_range_velocities[4]]),(w_velocity_nodes[1][x_pipe_range_velocities[1]]),(w_velocity_nodes[1][x_pipe_range_velocities[2] - 1]) , (w_velocity_nodes[1][x_pipe_range_velocities[3] + 1])], label = "pipe side walls", color = :deeppink2, linewidth = 1)
Colorbar(fig[2, 2], hm_T, label="C")
xS, yS, zS = nodes(S_t[1])
ax_S = Axis(fig[2, 3]; title="salinity[ppt]", axis_kwargs...)
hm_S = heatmap!(ax_S, xS, zS, Sₙ; colorrange=S_range, colormap=:haline) #note that this is still using old grid from T, S, initial, may need to recompute x and z using specific nodes 
xlims!(ax_S,(x_center - pipe_radius - pipe_wall_thickness - (1.5 * pipe_radius)), (x_center + pipe_radius + pipe_wall_thickness + (1.5 * pipe_radius)))
ylims!(ax_S, (-pipe_bottom_depth - (8 * pipe_radius)), (-pipe_top_depth + (8 * pipe_radius))) 
vlines!([(w_velocity_nodes[1][x_pipe_range_velocities[4]]),(w_velocity_nodes[1][x_pipe_range_velocities[1]]),(w_velocity_nodes[1][x_pipe_range_velocities[2] - 1]) , (w_velocity_nodes[1][x_pipe_range_velocities[3] + 1])], label = "pipe side walls", color = :deeppink2, linewidth = 1)
hideydecorations!(ax_S, grid = false)
Colorbar(fig[2, 4], hm_S, label="ppt")
xρ, yρ, zρ = nodes(ρ_t[1])
ax_ρ = Axis(fig[2, 5]; title="potential density[kg/m^3]", axis_kwargs...)
hm_ρ = heatmap!(ax_ρ, xρ, zρ, ρₙ; colorrange=ρ_range, colormap=Reverse(:viridis)) #note that this is still using old grid from T, S, initial, may need to recompute x and z using specific nodes 
xlims!(ax_ρ,(x_center - pipe_radius - pipe_wall_thickness - (1.5 * pipe_radius)), (x_center + pipe_radius + pipe_wall_thickness + (1.5 * pipe_radius)))
ylims!(ax_ρ, (-pipe_bottom_depth - (8 * pipe_radius)), (-pipe_top_depth + (8 * pipe_radius))) 
vlines!([(w_velocity_nodes[1][x_pipe_range_velocities[4]]),(w_velocity_nodes[1][x_pipe_range_velocities[1]]),(w_velocity_nodes[1][x_pipe_range_velocities[2] - 1]) , (w_velocity_nodes[1][x_pipe_range_velocities[3] + 1])], label = "pipe side walls", color = :deeppink2, linewidth = 1)
hideydecorations!(ax_ρ, grid = false)
Colorbar(fig[2, 6], hm_ρ, label="kg/m^3")
fig
@info "Making properties animation from data"
frames = 1:length(times)
record(fig, joinpath(pathname,"ZOOMproperties.mp4"), frames, framerate=8) do i
    @info string("Plotting frame ", i, " of ", frames[end])
    n[] = i
end


#uutility functions again 
#this function gives you an x index closest to x spacing coordinate - TODO: can make a lot faster via a search that halves things 

#this function gets you the max and min in pipe for an array 
function getMaxAndMinInPipeArr(arr, x_range::UnitRange, z_range::UnitRange)
    myMax = maximum(arr[x_range, z_range])
    myMin = minimum(arr[x_range, z_range])
    return (myMin, myMax)
end
#extract readable fields 
# ν = adapt(Array, interior(viscous_field, :, 1, :))
# b = adapt(Array, interior(buoyancy_field, :, 1, :))
#get max and min of the two to plot on the same scale
#not quite max and min in pipe, cuts out very top and bottom 
ν_range = getMaxAndMinInPipe(num_Data_Points, ν_t, (x_pipe_range_tracer[2]:x_pipe_range_tracer[3]), (getZIndex(w_velocity_nodes, -pipe_bottom_depth + 0.01)):(getZIndex(w_velocity_nodes, -pipe_top_depth - 0.01)))
b_range = getMaxAndMinInPipe(num_Data_Points, b_t, (x_pipe_range_tracer[2]:x_pipe_range_tracer[3]), (getZIndex(w_velocity_nodes, -pipe_bottom_depth + 0.01)):(getZIndex(w_velocity_nodes, -pipe_top_depth - 0.01)))
bν_colorbar_range = (-max(abs(ν_range[1]), abs(ν_range[2]), abs(b_range[1]), abs(b_range[2])), max(abs(ν_range[1]), abs(ν_range[2]), abs(b_range[1]), abs(b_range[2]))) #this makes the two nearly unreadable  


# plot them
x, y, z = nodes(ν_t[1])
fig = Figure(size = (1200, 700))
title = @lift @sprintf("t = %1.2f minutes", round(times[$n] / minute, digits=2))
axis_kwargs = (xlabel="x (m)", ylabel="z (m)", width=300)
fig[0, 1:4] = Label(fig, title)
#viscous component
ax_ν = Axis(fig[1, 1]; title="viscous component", axis_kwargs...)
# ν = adapt(Array, interior(viscous_field, :, 1, :))
hm_ν = heatmap!(ax_ν, x, z, νₙ; colorrange = bν_colorbar_range, colormap=:balance)
xlims!(ax_ν,(x_center - pipe_radius - pipe_wall_thickness - (0.5 * pipe_radius)), (x_center + pipe_radius + pipe_wall_thickness + (0.5 * pipe_radius)))
ylims!(ax_ν, (-pipe_bottom_depth - (8 * pipe_radius)), (-pipe_top_depth + (8 * pipe_radius)))  #note that this is still using old grid from T, S, initial, may need to recompute x and z using specific nodes 
vlines!([(w_velocity_nodes[1][x_pipe_range_velocities[4]]),(w_velocity_nodes[1][x_pipe_range_velocities[1]]),(w_velocity_nodes[1][x_pipe_range_velocities[2] - 1]) , (w_velocity_nodes[1][x_pipe_range_velocities[3] + 1])], label = "pipe side walls", color = :deeppink2, linewidth = 1)
# Colorbar(fig[1,2], hm_ν, label="m^2/s")
#buoyancy component 
ax_b = Axis(fig[1, 2]; title="Buoyancy component", axis_kwargs...)
hm_b = heatmap!(ax_b, x, z, bₙ; colorrange = bν_colorbar_range, colormap=:balance)
xlims!(ax_b,(x_center - pipe_radius - pipe_wall_thickness - (0.5 * pipe_radius)), (x_center + pipe_radius + pipe_wall_thickness + (0.5 * pipe_radius)))
ylims!(ax_b, (-pipe_bottom_depth - (8 * pipe_radius)), (-pipe_top_depth + (8 * pipe_radius)))  #note that this is still using old grid from T, S, initial, may need to recompute x and z using specific nodes 
vlines!([(w_velocity_nodes[1][x_pipe_range_velocities[4]]),(w_velocity_nodes[1][x_pipe_range_velocities[1]]),(w_velocity_nodes[1][x_pipe_range_velocities[2] - 1]) , (w_velocity_nodes[1][x_pipe_range_velocities[3] + 1])], label = "pipe side walls", color = :deeppink2, linewidth = 1)
hideydecorations!(ax_b, grid = false)
# Colorbar(fig[1,4], hm_b, label="m^2/s")
#buoyancy + viscous componnent
ax_bν = Axis(fig[1, 3]; title="Buoyancy + Viscous component", axis_kwargs...)
hm_bν = heatmap!(ax_bν, x, z, bνₙ; colorrange = bν_colorbar_range, colormap=:balance)
xlims!(ax_bν,(x_center - pipe_radius - pipe_wall_thickness - (0.5 * pipe_radius)), (x_center + pipe_radius + pipe_wall_thickness + (0.5 * pipe_radius)))
ylims!(ax_bν, (-pipe_bottom_depth - (8 * pipe_radius)), (-pipe_top_depth + (8 * pipe_radius)))  #note that this is still using old grid from T, S, initial, may need to recompute x and z using specific nodes 
vlines!([(w_velocity_nodes[1][x_pipe_range_velocities[4]]),(w_velocity_nodes[1][x_pipe_range_velocities[1]]),(w_velocity_nodes[1][x_pipe_range_velocities[2] - 1]) , (w_velocity_nodes[1][x_pipe_range_velocities[3] + 1])], label = "pipe side walls", color = :deeppink2, linewidth = 1)
hideydecorations!(ax_bν, grid = false)
Colorbar(fig[1,4], hm_bν, label="m^2/s")
fig
@info "Making viscosity and buoyancy heatmap animation"
frames = 1:length(times)
record(fig, joinpath(pathname,"componentsHeatmap.mp4"), frames, framerate=8) do i
    @info string("Plotting frame ", i, " of ", frames[end])
    n[] = i
end


#plot cross sections 
#extract cross sections
# ν_tenth =  ν[(x_pipe_range_tracer[1] - 4):(x_pipe_range_tracer[4] + 4) , z_index_tenth_tracer]
# ν_quarter =  ν[(x_pipe_range_tracer[1] - 4):(x_pipe_range_tracer[4] + 4) , z_index_quarter_tracer]
# ν_half =  ν[(x_pipe_range_tracer[1] - 4):(x_pipe_range_tracer[4] + 4) , z_index_half_tracer]
# ν_three_quarter=  ν[(x_pipe_range_tracer[1] - 4):(x_pipe_range_tracer[4] + 4) , z_index_three_quarter_tracer]
# b_tenth =  b[(x_pipe_range_tracer[1] - 4):(x_pipe_range_tracer[4] + 4) , z_index_tenth_tracer]
# b_quarter =  b[(x_pipe_range_tracer[1] - 4):(x_pipe_range_tracer[4] + 4) , z_index_quarter_tracer]
# b_half =  b[(x_pipe_range_tracer[1] - 4):(x_pipe_range_tracer[4] + 4) , z_index_half_tracer]
# b_three_quarter= b[(x_pipe_range_tracer[1] - 4):(x_pipe_range_tracer[4] + 4) , z_index_three_quarter_tracer]
# w_tenth =  interior(w_t[end], (x_pipe_range_velocities[1] - 4):(x_pipe_range_velocities[4] + 4) , 1, z_index_tenth_velocities)
# w_quarter =  interior(w_t[end], (x_pipe_range_velocities[1] - 4):(x_pipe_range_velocities[4] + 4) , 1, z_index_quarter_velocities)
# w_half=  interior(w_t[end], (x_pipe_range_velocities[1] - 4):(x_pipe_range_velocities[4] + 4) , 1, z_index_half_velocities)
# w_three_quarter =  interior(w_t[end], (x_pipe_range_velocities[1] - 4):(x_pipe_range_velocities[4] + 4) , 1, z_index_three_quarter_velocities)
νₙ_half =  @lift interior(ν_t[$n], (x_pipe_range_tracer[1] - 4):(x_pipe_range_tracer[4] + 4) , 1, z_index_half_tracer)
bₙ_half = @lift interior(b_t[$n], (x_pipe_range_tracer[1] - 4):(x_pipe_range_tracer[4] + 4) , 1, z_index_half_tracer)

νₙ_three_quarter =  @lift interior(ν_t[$n], (x_pipe_range_tracer[1] - 4):(x_pipe_range_tracer[4] + 4) , 1, z_index_three_quarter_tracer)
bₙ_three_quarter = @lift interior(b_t[$n], (x_pipe_range_tracer[1] - 4):(x_pipe_range_tracer[4] + 4) , 1, z_index_three_quarter_tracer)

fig = Figure(size = (1000, 600))
# myTicksBIG = -pipe_radius : pipe_radius/10: pipe_radius
# myTicksSMALL = -pipe_radius : pipe_radius/5: pipe_radius
kwargsBIG = (; xminorticks = IntervalsBetween(10), xminorticksvisible = true, xminorgridvisible = true)
kwargs =(; xminorticks = IntervalsBetween(5), xminorticksvisible = true, xminorgridvisible = true)
title = @lift @sprintf("t = %1.2f minutes", round(times[$n] / minute, digits=2))
fig[0, 1:4] = Label(fig, title)
cross_components_plot= Axis(fig[1,1:3], title = "w Momentum components @ 0.5 pipe", xlabel="x (m), centered at pipe center", ylabel = "acceleration m^2/s", width =  600; kwargsBIG...)
scatterlines!((x_plot_range_velocities .- x_center), νₙ_three_quarter, label = "viscous component", color = colors = :blue, markersize = 5)
scatterlines!((x_plot_range_velocities .- x_center), bₙ_three_quarter, label = "Buoyancy component", color = colors = :red, markersize = 5)
vlines!([(w_velocity_nodes[1][x_pipe_range_velocities[4]] - x_center),(w_velocity_nodes[1][x_pipe_range_velocities[1]] - x_center),(w_velocity_nodes[1][x_pipe_range_velocities[2] - 1] - x_center) , (w_velocity_nodes[1][x_pipe_range_velocities[3] + 1] - x_center)], label = "pipe side walls", color = :deeppink2, linewidth = 5)
ylims!(bν_colorbar_range[1], bν_colorbar_range[2])
cross_velocities_plot= Axis(fig[1,1:3], ylabel = "w velocity m/s", yaxisposition = :right)
hidespines!(cross_velocities_plot)
hidexdecorations!(cross_velocities_plot)
scatterlines!((x_plot_range_velocities .- x_center), wₙ_three_quarter , label = "w_velocity", color = colors = :black, markersize = 5)
ylims!(-max(abs(w_pipe_range[1]), abs(w_pipe_range[2])), max(abs(w_pipe_range[1]), abs(w_pipe_range[2])))
fig[1, 4] = Legend(fig, cross_components_plot, frame_visible = false)
fig
@info "Making viscosity and buoyancy plot animation"
frames = 1:length(times)
record(fig, joinpath(pathname,"componentsPlotThreeQuarter.mp4"), frames, framerate=8) do i
    @info string("Plotting frame ", i, " of ", frames[end])
    n[] = i
end



#FILTERED AVERAGES
averages_unfilt = zeros(length(times))
for j in 1 : length(times)
    averages_unfilt[j] = getMaskedAverageAtZ(pipeMask, w_t[j], (Center(), Center(), Face()), (-pipe_top_depth - 0.5*pipe_length), z_grid_spacing)
end
#filter
fs = 1/output_interval
# num_oscillation = 6
# fc = (1/num_oscillation) * (1/oscillation_period)
fc = 1/30minute
averages_filt = filtfilt(digitalfilter(Lowpass(fc, fs = fs), Butterworth(5)), averages_unfilt)
#plot
fig = Figure(size = (1000, 600))
title = "Velocity and Discharge"
fig[0, :] = Label(fig, title)
kwargs= (; yminorticks = IntervalsBetween(10), yminorticksvisible = true, yminorgridvisible = true)
velocities_plot= Axis(fig[1,1], title = "Pipe Velocities Averaged", xlabel="time(hrs)", ylabel = "velocities (m/s)", width =  700)
scatterlines!((times/hour), averages_unfilt, label = "velocity average min", color = :blue, markersize = 5)
scatterlines!((times/hour), averages_filt, label = "filtered velocity", color = :red, markersize = 0)
xlims!(0, times[end]/hour)
fig[1, 2] = Legend(fig, velocities_plot, frame_visible = false)
volume_flux_plot= Axis(fig[2, 1], title = "Discharge from filtered velocity & 3d conversion", xlabel="time(hrs)", ylabel = "discharge m^3/s", width =  700; kwargs...)
volume_flux = averages_filt .* (π*(pipe_radius)^2)
scatterlines!((times/hour), volume_flux, label = "discharge", color = :green, markersize = 0)
xlims!(0, times[end]/hour)
fig[2, 2] = Legend(fig, volume_flux_plot, frame_visible = false)
fig
save(joinpath(pathname,"Filtered Velocity and Discharge.png"), fig)
@info "Finished plotting average values vs time chart"


#plot lines showing velocity at different heights
#plot bar showing average flux through pipe













#the stuff below just plots one moment in time 
#GOAL: plot vorticity 
#some stuff for calculating the currently
# compute!(ν)
# compute!(b)
# νb_w = Field(ν + b)
# compute!(νb_w) #buoyancy and viscous field for w momentum

# u_viscousArgs = (model.closure, model.clock, model.velocities.u)
# u_viscous_component_kernel_op = KernelFunctionOperation{Center, Center, Center}(ViscousComponent, model.grid, u_viscousArgs...)
# ν_u = Field(u_viscous_component_kernel_op)
# compute!(ν_u) #vicous momentum equation in X

# ω = Field(-∂x(νb_w) + ∂z(ν_u)) ##angular 
# compute!(ω)


# x, y, z = nodes(ν_t[1])
# fig = Figure(size = (1200, 700))
# title = @lift @sprintf("t = %1.2f minutes", round(times[$n] / minute, digits=2))
# axis_kwargs = (xlabel="x (m)", ylabel="z (m)", width=300)
# fig[0, 1:4] = Label(fig, title)
# #w momentum 
# ax_ν = Axis(fig[1, 1]; title="w", axis_kwargs...)
# νb_w = adapt(Array, interior(νb_w, :, 1, :))
# hm_ν = heatmap!(ax_ν, x, z, νb_w; colorrange = bν_colorbar_range, colormap=:balance)
# xlims!(ax_ν,(x_center - pipe_radius - pipe_wall_thickness - (0.5 * pipe_radius)), (x_center + pipe_radius + pipe_wall_thickness + (0.5 * pipe_radius)))
# ylims!(ax_ν, (-pipe_bottom_depth - (8 * pipe_radius)), (-pipe_top_depth + (8 * pipe_radius)))  #note that this is still using old grid from T, S, initial, may need to recompute x and z using specific nodes 
# vlines!([(w_velocity_nodes[1][x_pipe_range_velocities[4]]),(w_velocity_nodes[1][x_pipe_range_velocities[1]]),(w_velocity_nodes[1][x_pipe_range_velocities[2] - 1]) , (w_velocity_nodes[1][x_pipe_range_velocities[3] + 1])], label = "pipe side walls", color = :deeppink2, linewidth = 1)
# # Colorbar(fig[1,2], hm_ν, label="m^2/s")
# #u momentum
# ax_b = Axis(fig[1,2]; title="u", axis_kwargs...)
# ν_u = adapt(Array, interior(ν_u, :, 1, :))
# hm_b = heatmap!(ax_b, x, z, ν_u; colorrange = bν_colorbar_range, colormap=:balance)
# xlims!(ax_b,(x_center - pipe_radius - pipe_wall_thickness - (0.5 * pipe_radius)), (x_center + pipe_radius + pipe_wall_thickness + (0.5 * pipe_radius)))
# ylims!(ax_b, (-pipe_bottom_depth - (8 * pipe_radius)), (-pipe_top_depth + (8 * pipe_radius)))  #note that this is still using old grid from T, S, initial, may need to recompute x and z using specific nodes 
# vlines!([(w_velocity_nodes[1][x_pipe_range_velocities[4]]),(w_velocity_nodes[1][x_pipe_range_velocities[1]]),(w_velocity_nodes[1][x_pipe_range_velocities[2] - 1]) , (w_velocity_nodes[1][x_pipe_range_velocities[3] + 1])], label = "pipe side walls", color = :deeppink2, linewidth = 1)
# hideydecorations!(ax_b, grid = false)
# Colorbar(fig[1,3], hm_bν, label="m^2/s")
# # Colorbar(fig[1,4], hm_b, label="m^2/s")
# #buoyancy + viscous componnent
# ax_bν = Axis(fig[1, 4]; title="curl of momentum equations", axis_kwargs...)
# ω = adapt(Array, interior(ω, :, 1, : ))
# ω_colorbar_range = getMaxAndMinInPipeArr(ω, x_pipe_range_tracer[2]:x_pipe_range_tracer[3], getZIndex(tracer_nodes, -pipe_bottom_depth + 0.01):getZIndex(tracer_nodes, -pipe_top_depth - 0.01))
# hm_bν = heatmap!(ax_bν, x, z, ω; colorrange = bν_colorbar_range, colormap=:balance)
# xlims!(ax_bν,(x_center - pipe_radius - pipe_wall_thickness - (0.5 * pipe_radius)), (x_center + pipe_radius + pipe_wall_thickness + (0.5 * pipe_radius)))
# ylims!(ax_bν, (-pipe_bottom_depth - (8 * pipe_radius)), (-pipe_top_depth + (8 * pipe_radius)))  #note that this is still using old grid from T, S, initial, may need to recompute x and z using specific nodes 
# vlines!([(w_velocity_nodes[1][x_pipe_range_velocities[4]]),(w_velocity_nodes[1][x_pipe_range_velocities[1]]),(w_velocity_nodes[1][x_pipe_range_velocities[2] - 1]) , (w_velocity_nodes[1][x_pipe_range_velocities[3] + 1])], label = "pipe side walls", color = :deeppink2, linewidth = 1)
# hideydecorations!(ax_bν, grid = false)
# Colorbar(fig[1,5], hm_bν, label="m^2/s")
# fig





















