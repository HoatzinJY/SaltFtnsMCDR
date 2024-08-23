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

const day = 86400;
const hour = 3600;
const minute = 60;

#search #SIMILITUDE to change things for similutde theories 
#IMPORTANT, IF WRITE DISCRETE FORCER HERE, NEED TO BE CAREFUL ABOUT CUARRAY VS ARRAY AND ADAPT THE CUARRAY OVER
"""NAME"""
trial_name = "SIMILITUDE to big"
mkdir(joinpath("Trials", (trial_name)))
pathname = joinpath("Trials", (trial_name))
filename = joinpath(pathname, "data")
#next run wih min 0.003, and then 2x
#then with 0.1, and 2x

#for full size model
#try 30m x 300m

"""COMPUTER parameters"""
const GPU_memory = 12

"""SIMULATION RUN INFORMATION"""
simulation_duration = 5hour
run_duration = 4minute
output_interval = 0.1
"""DOMAIN SIZE & SETUP"""
const domain_x = 0.07;#SIMILITUDE
const domain_z = 3.1; #SIMILITUDE
# const domain_z = 220; #BIG
const x_center = domain_x/2;
#max_grid_spacing = 0.02; #TODO: figure out what this needs to be set to 

"""PIPE SIZE AND SETUP"""
struct PipeWallData
    thermal_diffusivity :: Float64
    thickness :: Float64
end
const pipe_radius = 0.0015 #SIMILITUDE
const pipe_length = 2.7 #SIMILITUDE
const pipe_top_depth = 0.2 #SIMILITUDE
# const pipe_length = 200 #BIG
# const pipe_top_depth = 10 #BIG
const pipe_wall_thickness_intended = 0.0005 #Similitude, 0.005 orig
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
const sw_viscosity_molecular =  1.05e-6 * sqrt(5.4e-6)#SIMILITUDE
const sw_T_diffusivity_molecular = 1.46e-7 * sqrt(5.4e-6)#SIMILITUDE
#const sw_T_diffusivity_molecular = 1e-5 # as per the experimental data in papers Zhang 2004
const sw_S_diffusivity_molecular = 1.3e-9 
const sw_diffusivity_data = SeawaterDiffusivityData(sw_viscosity_molecular, sw_T_diffusivity_molecular, sw_S_diffusivity_molecular)


"""SET BACKGROUND TEMPERATURE AND SALINITY GRADIENTS"""
const T_top = 21.67;
const T_bot = 11.86;
const S_bot = 34.18;
const S_top = 35.22;
# const S_bot = 34.18;
# const S_top = 36.601543;
const delta_z = 200; 
const delta_z_multiplier = 1/2; #SIMILITUDE

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
roundUp(num::Float64, base) = ceil(Int, num / base) * base
findNearest(A::AbstractArray, x) = findmin(abs(A-x))
function getXYZ(i, j, k, grid, locs)
    myNodes = nodes(grid, locs, reshape = true)
    return[myNodes[1][i], myNodes[2][j], myNodes[3][k]]
end
function getMaxAndMin(numPoints, dataSeries)
    myMax = maximum(interior(dataSeries[1], :, 1, :))
    myMin = minimum(interior(dataSeries[1], :, 1, :))
    for i in 2:numPoints
        myMax = max(maximum(interior(dataSeries[i], :, 1, :)),myMax )
        myMin = min(minimum(interior(dataSeries[i], :, 1, :)), myMin)
    end
    return (myMin, myMax)
end
function getMaxAndMinArr(arr)
    myMax = maximum(arr)
    myMin = minimum(arr)
    return (myMin, myMax)
end
function getMaxAndMinEnd(dataSeries)
    myMax = maximum(interior(dataSeries[end], : , 1, :))
    myMin = minimum(interior(dataSeries[end], : , 1, :))
    return (myMin, myMax)
end
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
function getBackgroundDensity(z, Tfunc::Function, Sfunc::Function)#takes a positive depth
    eos = TEOS10EquationOfState() 
    S_abs = gsw_sa_from_sp(Sfunc(pipe_radius + 0.0001,z), gsw_p_from_z(z, 30), 31, -30) #random lat and long in ocean
    T_consv = gsw_ct_from_t(S_abs, Tfunc(pipe_radius + 0.0001,z), gsw_p_from_z(z, 30)) #set to 30 degrees north
    return TEOS10.ρ(T_consv, S_abs, z, eos) #note that this uses insitu s and T instead of conservative and absolute, which is what the function calls for 
end
function checkMemory(capacity, x_cells, z_cells)
    numCells = x_cells * z_cells
    possibleCells = (capacity/32) * 100000000
    if (numCells > possibleCells)
        throw("Too many cells for gpu memory")
        return numCells
    else 
        return numCells
    end
end     # uses the fact that 32gb is approx 100 million cells

"""CALCULATE DOMAIN DETAILS AND SAVE AS CONSTANTS"""
#find oscillation period to get relevant timescales & velocities
surrounding_density_gradient = (getBackgroundDensity(-pipe_top_depth - pipe_length, TWaterColumn, SWaterColumn) - getBackgroundDensity(-pipe_top_depth, TWaterColumn, SWaterColumn))/pipe_length #takes average for unaltered water column across the pipe (NOT DOMAIN!)
oscillation_angular_frequency =  sqrt((g/1000) * surrounding_density_gradient)
oscillation_period = 2π/oscillation_angular_frequency
@info @sprintf("N = %.3e rad/sec | Buoyancy Oscillation period: %.3f minutes",  oscillation_angular_frequency, oscillation_period/minute)
#grid spacing desired
max_grid_spacing = 0.0001 #SIMILITUDE
#max_grid_spacing = 0.95*((4 * sw_diffusivity_data.T * sw_diffusivity_data.ν)/(oscillation_angular_frequency^2))^(1/4)
@info @sprintf("Max Grid Spacing: %.3em", max_grid_spacing)
my_x_grid_spacing = max_grid_spacing; #max grid spacing is actually a bit of a misnomer, perhaps shoudl be better called min, its max in the sense that its the max resolution 
my_z_grid_spacing = max_grid_spacing;
#resolution
x_res = floor(Int, domain_x / my_x_grid_spacing);
z_res = floor(Int, domain_z / my_z_grid_spacing);
#grid spacing adjusted to be divisible in domain_grid
const x_grid_spacing = domain_x/x_res;
const z_grid_spacing = domain_z/z_res;
#call error if it is too large for model
checkMemory(GPU_memory, x_res, z_res)
#pipe wall thickness, as set by the model
const pipe_wall_thickness = roundUp(pipe_wall_thickness_intended, x_grid_spacing)
@info @sprintf("Pipe walls are %.3e meters thick", pipe_wall_thickness)
@info @sprintf("X spacings: %.3e meters | Z spacings: %.3e meters", x_grid_spacing, z_grid_spacing)
@info @sprintf("X resolution: %.3e | Z resolution: %.3e ", x_res, z_res)
#time stepping 
#location with equivalent density
#TODO: potentially fix that height displaced with an estimation of where the density is equal 
#v_max_predicted = oscillation_angular_frequency * 1
#0.03 stable
v_max_predicted = 0.001 #SIMILITUDE this is just based on old simulations, shoudl perhaps change this. #0.07 for 20m long, 20cm diameter pipe #0.001 for 1m long 0.02R pipe, 0.01 seems good for small lab scale
min_time_step_predicted = (min(x_grid_spacing, z_grid_spacing)*CFL)/v_max_predicted
max_time_step_allowed = 2 * min_time_step_predicted
#damping rate for forcers  
min_relaxation_timescale = 1.1 * max_time_step_allowed
max_relaxation_rate = 1/min_relaxation_timescale

"""MASKING FUNCTIONS"""
function noMask(x, z)
    return 1
end
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
function isPipeWall(x, z)
    left_wall_range = (x_center - pipe_radius - pipe_wall_thickness) .. (x_center - pipe_radius)
    right_wall_range = (x_center + pipe_radius) .. (x_center + pipe_radius + pipe_wall_thickness)
    vert_range = -(pipe_top_depth + pipe_length) .. -(pipe_top_depth)
    if ((in(x, left_wall_range) || in(x, right_wall_range)) && in(z, vert_range))
        return true
    elseif (isStopper(x, z))
        return true
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
function pipeWallVelocityMask(x, z)
    if (!isInsidePipe(x, z) && ((-pipe_bottom_depth) < z < (-pipe_top_depth)) && !isSidePipe(x, z))
        return 1
    elseif (isInsidePipe(x, z) && isPipeWall(x, z)) #this is for the stopper that is for flow sepparation 
        return 1
    else
        return 0
    end
end

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
function tracerRelaxationMaskDomainTwo(x, y, z)
    #if not pipe wall, pipe, surface, or two sides
    unmasked_ratio = 0.4
    if (!isPipeWall(x, y, z) && !isInsidePipe(x, y, z) && (z < (-pipe_top_depth)) && !isSidePipe(x, y, z))
        return 1
    elseif (z > (-pipe_top_depth + unmasked_ratio * pipe_top_depth))
        return (0.5/((1 - unmasked_ratio)*pipe_top_depth))*z + 0.5
    else
        return 0
    end
end
function velocityRelaxationMaskDomainOne(x, y, z)
    if (!isPipeWall(x, y, z) && !isInsidePipe(x, y, z) && ((-pipe_bottom_depth) < z < (-pipe_top_depth)) && !isSidePipe(x, y, z))
        return 1
    else
        return 0
    end
end 
function tracerRelaxationMaskDomainThree(x, y, z)
    unmasked_ratio = 0.4
    if (min(x, domain_x - x) < (domain_z/2) - pipe_radius - pipe_wall_thickness - (pipe_length/10))
        return 1
    elseif (z > (-pipe_top_depth + unmasked_ratio * pipe_top_depth))
        return (0.5/((1 - unmasked_ratio)*pipe_top_depth))*z + 0.5
    elseif (z < (-pipe_bottom_depth - unmasked_ratio * abs(domain_z - pipe_bottom_depth)))
        return (0.5/((1 - unmasked_ratio)*(domain_z - pipe_bottom_depth)))*(-z - pipe_bottom_depth - ((1 - unmasked_ratio)*(domain_z - pipe_bottom_depth)))
    else
        return 0
    end
end
isInsidePipe(x, y, z) = isInsidePipe(x, z)
pipeMask(x,y,z) = pipeMask(x, z)
isSidePipe(x, y, z) = isSidePipe(x, z)
sidePipeMask(x, y, z) = sidePipeMask(x, z)
isPipeWall(x, y, z) = isPipeWall(x, z)
pipeWallMask(x, y, z) = pipeWallMask(x, z)
isSidePipeWall(x, y, z) = isSidePipeWall(x, z)
sidePipeWallMask(x, y, z) = sidePipeWallMask(x, z)
tracerRelaxationMaskDomainTwo(x, z) = tracerRelaxationMaskDomainTwo(x, 0, z)
velocityRelaxationMaskDomainOne(x, z) = velocityRelaxationMaskDomainOne(x, 0, z)
tracerRelaxationMaskDomainThree(x, z) = tracerRelaxationMaskDomainThree(x, 0, z)
pipeWallVelocityMask(x, y, z) = pipeWallVelocitymask(x, z)

"""INITIAL CONDITIONS"""
#this function for some reason does not work, calculates a gradient that is too large 
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
    return ((ρ₀ - ρ_wc_avg))/(α_avg * 1000 * p_length)
end 

# const pipeTGrad = neutralDensityPipeTempGradient(pipe_length, pipe_top_depth, max_grid_spacing, SWaterColumn, TWaterColumn)
const pipeTGrad = 0.0206#similitude 0.0418 for 1m long 

 #const pipeTGrad = 0.02

function T_init(x, y, z)
    if (isInsidePipe(x, z) || isPipeWall(x, z))
        return TWaterColumn(-pipe_bottom_depth) + (pipeTGrad * (z + pipe_bottom_depth))
    elseif (isSidePipe(x, z))
        return TWaterColumn(-pipe_top_depth) + (pipeTGrad * (z + pipe_top_depth))
    else
        return TWaterColumn(z)
    end
end


#fully equillibriated 
# function T_init(x, y, z)
#     return TWaterColumn(z)
# end

function S_init(x, y, z)
    if (isInsidePipe(x, z) || isPipeWall(x, z))
        return SWaterColumn(-pipe_bottom_depth)
    elseif (isSidePipe(x, z))
        return SWaterColumn(-pipe_top_depth)
    else
        return SWaterColumn(z)
    end
end
function A_init(x, y, z)
    if(z < -pipe_bottom_depth)
        return 1
    else
        return 0
    end
end 
function w_init(x, z)
    return 0;
end
T_init(x, z) = T_init(x, 0, z)
S_init(x, z) = S_init(x, 0, z)
A_init(x, z) = A_init(x, 0, z)

#diagnostic function to check buoyancy differences
function getBackgroundDensityPipe(z, Tfunc::Function, Sfunc::Function)#takes a positive depth
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
    wc_density_tot += getBackgroundDensity(z, T_init, S_init)
    pipe_density_tot += getBackgroundDensityPipe(z, T_init, S_init)
    equalized_pipe_density_tot += getBackgroundDensityPipe(z, TWaterColumn, S_init) #as if fully equalized, for similitude 
    z += z_grid_spacing
    count += 1
end

@info @sprintf("wc density - pipe_density = %.2e kg/m^3", (wc_density_tot/count) - (pipe_density_tot/count))
@info @sprintf("@ eq, wc density - pipe_eq_density = %.3ek kg/m^3", (wc_density_tot/count) - (equalized_pipe_density_tot/count))

"""MODEL BOUNDARY CONDITIONS SETUP"""
# get gradients for boundary conditons
initial_T_top_gradient = (TWaterColumn(0) - TWaterColumn(0 - z_grid_spacing))/z_grid_spacing
initial_T_bottom_gradient = (TWaterColumn(domain_z + z_grid_spacing) - TWaterColumn(0domain_z))/z_grid_spacing
initial_S_top_gradient = (SWaterColumn(0, 0 - z_grid_spacing) - SWaterColumn(0,0))/z_grid_spacing
initial_S_bottom_gradient = (SWaterColumn(0, domain_z) - SWaterColumn(0, domain_z + z_grid_spacing))/z_grid_spacing

"""FORCING FUNCTIONS"""
# #none currently used 
TWaterColumnTarget(x, z, t) = TWaterColumn(z)
SWaterColumnTarget(x, z, t) = SWaterColumn(z)

"""MODEL DIFFUSIVITY SETUP"""
const diffusivity_data = (seawater = sw_diffusivity_data, pipe = pipe_data)
const pipeWallThickness = (actual = pipe_wall_thickness, intended = pipe_wall_thickness_intended)
function tempDiffusivities(x, y, z, diffusivities::NamedTuple, wallThickness::NamedTuple, wall_indicator::String)
    return diffusivities[:seawater].T
    # # if (isPipeWall(x, z) || wall_indicator == "WALL")
    # #     return diffusivities[:pipe].thermal_diffusivity * (wallThickness[:actual]/wallThickness[:intended]) #edit to account for wall thickness
    # # else 
    #     #try scaling 
    #     max_velocity_predicted = 0.03
    #     return max_velocity_predicted*max_grid_spacing*1.1
    #     #return diffusivities[:seawater].T
    # # end
end
#crashed with zeroing out salt diffusivity and viscosity
function saltDiffusivities(x, y, z, diffusivities::NamedTuple)
    if(velocityRelaxationMaskDomainOne(x, y, z) == 1)
        return diffusivities[:seawater].S
    elseif (isPipeWall(x, z) || isSidePipeWall(x, z))
        return 0
    else
        return diffusivities[:seawater].S
    end
end
function myViscosity(x, y, z, diffusivities::NamedTuple)
    if(velocityRelaxationMaskDomainOne(x, y, z) == 1)
        return diffusivities[:seawater].ν
    elseif (isPipeWall(x, z) || isSidePipeWall(x, z))
        return 0
    else 
        return diffusivities[:seawater].ν
    end
end
tempDiffusivities(x, z, t) = tempDiffusivities(x, 0, z, diffusivity_data, pipeWallThickness, "")
saltDiffusivities(x ,z ,t ) = saltDiffusivities(x, 0, z, diffusivity_data)
myViscosity(x ,z ,t ) = myViscosity(x, 0, z, diffusivity_data)

"""MODEL SETUP"""
domain_grid = RectilinearGrid(GPU(), Float64; size=(x_res, z_res), x=(0, domain_x), z=(-domain_z, 0), topology=(Periodic, Flat, Bounded))

clock = Clock{eltype(domain_grid)}(time=0)

timestepper = :QuasiAdamsBashforth2

advection = WENO()

eos = TEOS10EquationOfState()
myBuoyancy = SeawaterBuoyancy(equation_of_state=eos)


tracers = (:T, :S, :A)
closure = ScalarDiffusivity(ν=myViscosity, κ=(T=tempDiffusivities, S=saltDiffusivities, A = 0))
#closure = ScalarDiffusivity(ν = 1, κ=1)

#domain is periodic in the east-west
T_bcs = FieldBoundaryConditions(top = GradientBoundaryCondition(initial_T_top_gradient), bottom = GradientBoundaryCondition(initial_T_bottom_gradient))
S_bcs = FieldBoundaryConditions(top = GradientBoundaryCondition(initial_S_top_gradient), bottom = GradientBoundaryCondition(initial_S_bottom_gradient))
boundary_conditions = (T = T_bcs, S = S_bcs)

#TODO: combine these when set 
#uw_domain_forcer = Relaxation(rate = max_relaxation_rate, mask = velocityRelaxationMaskDomainOne)
uw_pipe_wall_forcer = Relaxation(rate = max_relaxation_rate, mask = pipeWallVelocityMask)
T_domain_forcer = Relaxation(rate = max_relaxation_rate, mask = tracerRelaxationMaskDomainTwo, target = TWaterColumnTarget)
S_domain_forcer = Relaxation(rate = max_relaxation_rate, mask = tracerRelaxationMaskDomainTwo, target = SWaterColumnTarget)
#S_side_forcer = Relaxation(rate = max_relaxation_rate, mask = sidePipeMask, target = SWaterColumn(-pipe_top_depth))
forcing = (u = (uw_pipe_wall_forcer),  w = (uw_pipe_wall_forcer), T = (T_domain_forcer), S = (S_domain_forcer))

model = NonhydrostaticModel(; grid=domain_grid, clock, advection, buoyancy = myBuoyancy, tracers, timestepper, closure, forcing, boundary_conditions)
# model.output_writers[:checkpointer] = Checkpointer(model; frequency = 150000, dir = pathname, prefix = "checkpoint", force = false, verbose = false, properties = [:architecture, :boundary_conditions, :grid, :clock, :buoyancy, :closure, :velocities, :tracers, :timestepper])

# restore_iteration  = 150000
# checkpoint_filename = joinpath(pathname, "checkpoint"*string(restore_iteration)) #TODO: check this 
# model_kwargs = Dict("closure" => closure, "forcing" => forcing)
# model = restore_from_checkpoint(checkpoint_filename; model_kwargs...)

density_operation = seawater_density(model; geopotential_height)
@info "model made"


"""SET UP INITIAL CONDITIONS"""
set!(model, T= T_init, S=S_init, w = w_init, A = A_init)#ASK  - why set twice ??
set!(model, T= T_init, S=S_init, w = w_init, A = A_init)
@info "initial conditions set"


"""SETTING UP SIMULATION"""
#finding various time scales 
min_grid_spacing = min(minimum_xspacing(model.grid), minimum_zspacing(model.grid))
initial_travel_velocity_fake = 0.025 # a number just to set an initial time step, i set it to be around max for safety
initial_advection_time_scale = min_grid_spacing/initial_travel_velocity_fake
diffusion_time_scale = (min_grid_spacing^2)/model.closure.κ.T(0, 0, 0, diffusivity_data, pipeWallThickness, "WALL") #TRACER_MIN, set to tracer with biggest kappa, to use when using function based diffusivity
# diffusion_time_scale = (min_grid_spacing^2)/model.closure.κ.T
viscous_time_scale = (min_grid_spacing^2)/model.closure.ν(0, 0, 0, diffusivity_data)
#initial_time_step = 0.5*min(0.2 * min(diffusion_time_scale, oscillation_period, viscous_time_scale, initial_advection_time_scale), max_time_step_allowed)
initial_time_step = 0.0048 #0.0048 for 1m 0.02R
"""IMPORTANT, diffusion cfl does not work with functional kappas, need to manually set max step"""
new_max_time_step = min(0.2 * diffusion_time_scale, 0.2 * viscous_time_scale, max_time_step_allowed) #TRACER_MIN, uses a cfl of 0.2, put in diffusion time scale of tracer with biggest kappa, or viscosity
#set up simulation & timewizard
simulation = Simulation(model, Δt=initial_time_step, stop_time=simulation_duration, wall_time_limit=run_duration) # make initial delta t bigger



#various callbacks
#timewizard
timeWizard = TimeStepWizard(cfl=CFL, diffusive_cfl = CFL, max_Δt = new_max_time_step, max_change = 1.01, min_change = 0.5) 
simulation.callbacks[:timeWizard] = Callback(timeWizard, IterationInterval(4))
#progress Message
progress_message(sim) = @printf("Iteration: %04d, time: %s, Δt: %s, wall time: %s\n", iteration(sim), prettytime(sim), prettytime(sim.Δt), prettytime(sim.run_wall_time))
add_callback!(simulation, progress_message, IterationInterval(50))

#set up output_writers for large fields
w = model.velocities.w;
u = model.velocities.u;
T = model.tracers.T;
S = model.tracers.S;
ζ = Field(-∂x(w) + ∂z(u)) #vorticity in y 
ρ = Field(density_operation)
A = model.tracers.A
function getXIndex(fieldNodes, xCoord)
    for i in 1 : length(fieldNodes[1])
        if(abs(fieldNodes[1][i] - xCoord) <= 1.01(x_grid_spacing/2)) #1.10 adjusts for some rounding issue with centers 
            return i
        end
    end
    throw("index not found for given z value")
end
@inline function ViscousComponent(i, j, k, grid, closure, clock, velocities)
    #return _νᶜᶜᶜ(i, j, k, grid, model.closure, model.diffusivity_fields, model.clock)  * ∇²ᶜᶜᶜ(i, j, k, grid, model.velocities.w)
    return _νᶜᶜᶜ(i, j, k, grid, closure, nothing, clock)  * ∇²ᶜᶜᶜ(i, j, k, grid, velocities)
end
w_viscousArgs = (model.closure, model.clock, model.velocities.w)
w_viscous_component_kernel_op = KernelFunctionOperation{Center, Center, Center}(ViscousComponent, model.grid, w_viscousArgs...)
ν= Field(w_viscous_component_kernel_op)
#viscous_field = compute!(ν_component)
#this computes the buoyancy relative to the point stationary outside --> not right, but how?
@inline function wcBuoyancy(i, j, k, grid, b::SeawaterBuoyancy, tracer_fields)
    T, S = get_temperature_and_salinity(b, tracer_fields)
    return (- (b.gravitational_acceleration * ρ′(i, j, k, grid, b.equation_of_state, T, S)
    / b.equation_of_state.reference_density))
end
#this computes the buoyancy acceleration relative to a reference value
@inline function BuoyancyComponent(i, j, k, grid, b::SeawaterBuoyancy, tracer_fields, wc_index)
    T, S = get_temperature_and_salinity(b, tracer_fields)
    return (- (b.gravitational_acceleration * ρ′(i, j, k, grid, b.equation_of_state, T, S)
    / b.equation_of_state.reference_density)) - wcBuoyancy(wc_index, j, k, grid, b, tracer_fields)
end
tracer_nodes = nodes(model.grid, (Center(), Center(), Center()), reshape = true)
wc_index = getXIndex(tracer_nodes, pipe_radius + (3 *x_grid_spacing))
buoyancyArgs = (myBuoyancy, model.tracers, wc_index) #for some reason model.buoyancy.equation_of_state does not return an eos
#buoyancy_component_kernel_op = KernelFunctionOperation{Center, Center, Center}(buoyancy_perturbationᶜᶜᶜ, model.grid, buoyancyArgs...) #uses built in
buoyancy_component_kernel_op = KernelFunctionOperation{Center, Center, Center}(BuoyancyComponent, model.grid, buoyancyArgs...)
b = Field(buoyancy_component_kernel_op)
#buoyancy_field = compute!(b_component)






simulation.output_writers[:outputs] = JLD2OutputWriter(model, (; u, w, T, S, ζ, ρ, A, ν, b); filename, schedule=TimeInterval(output_interval), overwrite_existing=true) #can also set to TimeInterval
#simulation.output_writers[:outputs] = JLD2OutputWriter(model, (; u, w, T, S, ζ, ρ, A); filename, schedule=TimeInterval(output_interval), overwrite_existing=true) #can also set to TimeInterval


run!(simulation, pickup = false)
# simulation.wall_time_limit +=10hour
# run!(simulation)



#visualize simulation
#visualize simulation
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





















