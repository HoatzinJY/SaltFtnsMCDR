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

"""TODO: NEW GRID SPACING AND MAX TIME STEP WITHOUT WORRYING ABOUT OSCILLATION AND INSTEAD ABOUT BOUNDARY LAYER, 
PERHAPS SHOULD KEEP CURRENT BECUASE THERE IS OSCILLATION AT THE END
SHOUDL REALLY THINK ABOUT THE TIMESTEP PROBLEM, CANT BE RELAXING TOO SLOW, MAYBE ADD A CALLBACK ON RELAXATION TIME
TO BE A FUNCTION OF THE TIMESTEP?"""

#IF ADD IN MORE TRACERS OR HAVE WALL SPECIFIC ONES, LOOK FOR LINES LABELED TRACER_MIN, need to edit, add function that takes min/max of things
#to swap between bounded and periodic, look for LAB_SET
#NEXT RUN: NO WALLS TRIAL
# EVEN WITH CORRECT PERTURBATION IT DIVERGES
# realized: STABLE WITH 10x steeper gradient but not this, also unstable once a sponge layer is thrown in 
#thought: looking at traceradvection term, look at how time step wizard does normal cfl 
#thought: add nutrients, carbon
#think --> look at diffusion, and relaxation term, how to get target from field 
#think --> look at how to set variable kappa
#TODO: figure out field argument and how to manipulate 
#when reading file, maybe use discrete?
#= 
- add in ability to take in external field
- make wall forcings accurate to friction
- make wall forcings accurate to diffusion  and block salt(think about equation, change in tracer with c, equate that to flux of heat through pipe walls)
- make diffusion and viscous terms accurate to mixing 
- accurately model infinity reservoir
- figure out issue with temperature --> potentially graph the average field? 
# TODO: to resolve above problem, maybe average temperature over certain grid cells
to do above, can use 
Oceananigans.AbstractOperations.Average, dimensions is three tuple,
=#
#TODO: problem: wall
#TODO: organizize file, declare variables as global and make function file, change init to profile 
#constants and operations
#TODO: figure out why it seems to be diverging, timestep too large, or, diverges over approx 50min
#linked with above, dampening?, should besomething to do with viscosity? 
#TODO: accept external temp files


"""CONSTANTS AND PHYSICAL PROPERTIES"""
const day = 86400;
const hour = 3600;
const minute = 60;
const g = 9.806; 
const CFL = 0.2;#gravitational acceleration, in m^2/s
mutable struct DiffusionParameter
    molecular::Float64
    isopycnal::Float64 #horizontal
    diapycnal::Float64 #vertical
end
mutable struct SeawaterDiffusivityData
    T :: DiffusionParameter
    S :: DiffusionParameter
end
#numbers from DPO, at 20C, low estimate for isopycnal
#TODO: maybe make values store a range for temp, and then search in range after interpolation
eddy_horizontal_diffusivity = 5e2
eddy_vertical_diffusivity = 1e-5
my_viscosity = DiffusionParameter(1.05e-6, 1e3, 1e-4)
T_diffusivity = DiffusionParameter(1.46e-7, eddy_horizontal_diffusivity, eddy_vertical_diffusivity )
S_diffusivity = DiffusionParameter(1.3E-9, eddy_horizontal_diffusivity, eddy_vertical_diffusivity )
seawater_diffusion_data = SeawaterDiffusivityData(T_diffusivity, S_diffusivity)
geopotential_height = 0; # sea surface height for potential density calculations
#useful auxilliary functions



"""COMPUTATIONAL FUNCTIONS"""
roundUp(num::Float64, base) = ceil(Int, num / base) * base
findNearest(A::AbstractArray, x) = findmin(abs(A-x)) # returns [nearest value, value index]
sigmoid(var, center, steepness) = 1/(1 + exp(-(var - center)/steepness))
xFlippedSigmoid(var, center, steepness) = 1/(1 + exp((var - center)/steepness))
function getMaxAndMin(numPoints, dataSeries)
    myMax = maximum(interior(dataSeries[1], :, 1, :))
    myMin = minimum(interior(dataSeries[1], :, 1, :))
    for i in 2:numPoints
        myMax = max(maximum(interior(dataSeries[i], :, 1, :)),myMax )
        myMin = min(minimum(interior(dataSeries[i], :, 1, :)), myMin)
    end
    return (myMin, myMax)
end
function getIndexFromTime(time, timeSeries)
    numTotFrames = length(timeSeries)
    if (time == -1)
        return numTotFrames
    end
    for i in 1:numTotFrames
        if (timeSeries[i] > time)
            return i
        end
    end
    @info "Run duration too short, animating until end"
    return numTotFrames
end 
#gets average over a masked area
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
#returns background density of the general water column at depth z (before perturbation)
function getBackgroundDensity(z, Tfunc::Function, Sfunc::Function)#takes a positive depth
    eos = TEOS10EquationOfState() 
    S_abs = gsw_sa_from_sp(Sfunc(0,z), gsw_p_from_z(z, 30), 31, -30) #random lat and long in ocean
    T_consv = gsw_ct_from_t(S_abs, Tfunc(0,z), gsw_p_from_z(z, 30)) #set to 30 degrees north
    return TEOS10.ρ(T_consv, S_abs, z, eos) #note that this uses insitu s and T instead of conservative and absolute, which is what the function calls for 
end
# """REAL SIZE ISH TEST"""
# scale = 1
# domain_x = 40scale # width, in meters
# domain_z = 300scale #height, in meters
# x_center = domain_x / 2;
# pipe_radius = 1scale;
# pipe_length = 200scale;
# pipe_top_depth = 50scale;
# pipe_bottom_depth = pipe_top_depth + pipe_length
# pipe_wall_thickness_intended = 0.02scale; #will be rounded up to nearest number of grid cells during grid setup
# #pumping parameters
# height_displaced = 40scale;
# initial_pipe_velocity = 0.001scale;



"""DOMAIN SIZE"""
scale = 1
domain_x = 15scale # width, in meters
domain_z = 18scale #height, in meters
x_center = domain_x / 2;



"""PIPE AND PUMPING PARAMETERS"""
#TODO: make these a struct 
pipe_radius = 0.5scale;
pipe_length = 10scale;
pipe_top_depth = 3scale;
pipe_bottom_depth = pipe_top_depth + pipe_length
pipe_wall_thickness_intended = 0.01scale; #will be rounded up to nearest number of grid cells during grid setup
#pumping parameters
height_displaced = 1scale;
initial_pipe_velocity = 0.001scale; #not yet used
#pipe wall properties 
mutable struct PipeWallData
    wall_thermal_diffusivity :: Float64 
    wall_thickness ::Float64
end
#current material: 90/10 copper nickel alloy https://super-metals.com/wp-content/uploads/2015/03/CuNI-90-10.pdf
wall_material_ρ = 8900#density, in kg/m^3
wall_material_cₚ = 376.812 #specific heat capacity under constant pressure, J/(kg * K)
wall_material_k = 50.208#thermal conductivity (W/(mK))
pipe_data = PipeWallData(wall_material_k/(wall_material_cₚ*wall_material_ρ), pipe_wall_thickness_intended)
@info @sprintf("Wall Material Diffusivity Unscaled: %.3e", pipe_data.wall_thermal_diffusivity)



"""MASKING FUNCTIONS"""
#TODO: these all use a pipe_wall_thickness --> should pass that in instead, since that is calculated during grid making
function isInsidePipe(x, z)
    if (x > (x_center - pipe_radius) && x < (x_center + pipe_radius) && z < -pipe_top_depth && z > -(pipe_top_depth + pipe_length))
        return true
    else
        return false
    end
end
isInsidePipe(x, y, z) = isInsidePipe(x, z)
function pipeMask(x, z)
    if (isInsidePipe(x, z))
        return 1
    else
        return 0
    end
end
pipeMask(x,y,z) = pipeMask(x, z)
#returns true or false for if coordinate is inside imaginary pipe walls, relies on center grid 
function isPipeWall(x, z)
    left_wall_range = (x_center - pipe_radius - pipe_wall_thickness) .. (x_center - pipe_radius)
    right_wall_range = (x_center + pipe_radius) .. (x_center + pipe_radius + pipe_wall_thickness)
    vert_range = -(pipe_top_depth + pipe_length) .. -(pipe_top_depth)
    if ((in(x, left_wall_range) || in(x, right_wall_range)) && in(z, vert_range))
        return true
    else
        return false
    end
end
isPipeWall(x, y, z) = isPipeWall(x, z)
function pipeWallMask(x,z)
    if (isPipeWall(x, z))
        return 1
    else
        return 0
    end
end
pipeWallMask(x, y, z) = pipeWallMask(x, z)
function noMask(x, z)
    return 1
end
#no mask
noMask(x, y, z) = noMask(x, z)
#border sponge layer 
function borderMask(x, z)
        if(min(x, domain_x - x) < x_grid_spacing || min(-z, domain_z + z) < z_grid_spacing)
            return 1
        else
            return 0
        end
end
borderMask(x, y, z) = borderMask(x, z)
#border mask on borders that are in contact with the water, 2 grid cells thick
function  waterBorderMask(x, z)
    if(min(x, domain_x - x) < 2(x_grid_spacing) || (domain_z + z) < 2(z_grid_spacing))
        return 1
    else
        return 0
    end
end
waterBorderMask(x, y, z) = waterBorderMask(x, z)
#to mimic prevent diffusion through pipe walls
#TODO: change around these constants and see
#TODO: figure out how toinline functions
#this function requires the grid to first be processes to give it actual pipe wall width, and for the pipe top depth to be given 
function tracerRelaxationMask(x, y, z)
    #no masking above top of pipe
    if(z > (-pipe_top_depth))
        return 0
    elseif(isInsidePipe(x, y, z) || isPipeWall(x, y, z))
        return 0
    elseif(z > (-pipe_top_depth - (pipe_length/2))) #top half of pipe 
        #sets of constants to play with, TODO: delete section once set, for efficiency 
        #z_perturb_base = -3.5 #uses z sign convention
        #z_steepness = 0.3 #lower numbers = more steep
        #TODO: calculate neccessary steepness to have mask of approximately 1 by the time 1 grid cell has happened
        #TODO: make these a function of grid size? 
        x_steepness_final = 5 #will be slightly larger than this by x steepness initial
        x_steepness_initial = 0.01 #the smaller this is, the steeper
        xz_steepness = 5 #how quickly x mask becomes gentler as you go down, smaller = quicker
        x_steepness = x_steepness_final*2 * sigmoid((-z),(pipe_top_depth),xz_steepness) - (x_steepness_final) + x_steepness_initial
        if(x > x_center)
            x_right_mask =  2 * xFlippedSigmoid(x, (x_center + pipe_radius + pipe_wall_thickness),x_steepness)
            #z_mask = @sigmoid((-z), (-z_perturb_base), z_steepness)
            return x_right_mask
        else 
            x_left_mask = 2 * sigmoid(x, (x_center - pipe_radius - pipe_wall_thickness), x_steepness)
            #z_mask = @sigmoid((-z), (-z_perturb_base), z_steepness)
            return x_left_mask
        end 
    else #for bottom half of pipe, return mirror image of top half 
        mirror_depth = -(pipe_top_depth + pipe_length/2) + ((-z) - (pipe_top_depth + pipe_length/2))
        return tracerRelaxationMask(x, y, mirror_depth)
    end
end
#an option that sets just the layer immediately to the side to be forced 
#function requires a grid size to already be set
function tracerRelaxationMaskStrip(x, y, z)
    if (x > x_center) #right side
        if (isPipeWall(x - x_grid_spacing, z))
            return 1
        else
            return 0
        end
    else #left side
        if(isPipeWall(x + x_grid_spacing, z))
            return 1
        else
            return 0
        end
    end
end
tracerRelaxationMask(x, z) = tracerRelaxationMask(x, 0, z)
tracerRelaxationMaskStrip(x, z) = tracerRelaxationMaskStrip(x, 0, z)



"""INITIAL CONDITIONS & WATER COLUMN CONDITIONS"""
#sets initial t and s. Assumes that z will be negative. 
#TODO:set up container functions for these to allow moving to another file 
#TODO: rename to TInit
T_top = 21.67;
T_bot = 11.86;
S_bot = 34.18;
S_top = 35.22;
delta_z = 200; 
#no perturbation
# function T_init(x, y, z)
#         return T_top - ((T_bot - T_top) / delta_z)z
# end
# function S_init(x, y, z)
#         return S_top - ((S_bot - S_top) / delta_z)z
# end
#perturbation
# function T_init(x, y, z)
#     #inside pipe
#     if (isInsidePipe(x, z))
#         return T_top - ((T_bot - T_top) / delta_z) * (z - height_displaced)
#         #outside pipe
#     else
#         return T_top - ((T_bot - T_top) / delta_z)z
#     end
# end
# function S_init(x, y, z)
#     if (isInsidePipe(x, z))
#         return S_top - ((S_bot - S_top) / delta_z) * (z - height_displaced)
#     else
#         return S_top - ((S_bot - S_top) / delta_z)z
#     end
# end
#TODO: test these initial conditions
#assuming having equalized temperature, steady state
function T_init(x, y, z)
    return T_top - ((T_bot - T_top) / delta_z)z
end
#this function does not need pipe, just the background profile, which
function S_init_b(x, y, z)
    return S_top - ((S_bot - S_top) / delta_z)z
end
function S_init(x, y, z)
    if (isInsidePipe(x, z) || isPipeWall(x, z))
        return S_top - ((S_bot - S_top) / delta_z)*(-pipe_bottom_depth)
    else
        return S_top - ((S_bot - S_top) / delta_z)z
    end
end
#assuming not having equalized temperatures, but salinity equal to bottom
# function T_init(x, y, z)
#     #inside pipe, currently slightly lagging by 0.2 m
#     if (isInsidePipe(x, z))
#         return T_top - ((T_bot - T_top) / delta_z) * (z - height_displaced)
#         #outside pipe
#     else
#         return T_top - ((T_bot - T_top) / delta_z)z
#     end
# end
# function S_init(x, y, z)
#     if (isInsidePipe(x, z))
#         return S_top - ((S_bot - S_top) / delta_z)*(-pipe_bottom_depth)
#     else
#         return S_top - ((S_bot - S_top) / delta_z)z
#     end
# end
T_init(x, z) = T_init(x, 0, z)
T_init_bc(y, z, t) = T_init(0, 0, z)
T_init_bc(z, t) = T_init(0, z)
T_init_target(x, z, t) = T_init(x, 0, z)
T_init_target(x, y, z, t) = T_init(x, 0, z)
S_init(x, z) = S_init(x, 0, z)
S_init_b(x, z) = S_init_b(x, 0, z)
S_init_bc(y, z, t) = S_init(0, 0, z)
S_init_bc(z, t) = S_init(0, z)
S_init_target(x, z, t) = S_init(x, 0, z)
S_init_target(x, y, z, t) = S_init(x, 0, z)
function w_init(x, y, z)
    if (isInsidePipe(x, z))
        #return initial_pipe_velocity;
        return 0
    else
        return 0
    end
end
w_init(x, z) = w_init(x, 0, z)
#TODO: get a word for "general water column"  
"""FORCING FUNCTIONS"""
# #TODO: test this 
function w_pump(x, y, z, t, w, p)
    pump_time = p.Δz / p.wₚ
    if(isInsidePipe(x, z) && t < pump_time)
        return p.f_rate * (p.wₚ - w)
    else
        return 0
    end
end
w_pump(x, z, t, w, p) = w_pump(x, 1, z, t, w, p)
# #TODO: test this
#this function forces the values in the wall to be equal to that immediately inside the pipe 
#takes in discrete form location i, j, k, a model field, and a relaxation rate \
#TODO: this function could potentially be really slow, maybe 
#TODO: test wall forcer --> tested manually, checked that it correctly identifies walls and gives correct forcing term
function inpenetrable_wall_forcer(x_center, fieldNodes, i, j, k, field, rate)
    # on left side of pipe
    if(fieldNodes[1][i] < x_center)
        return @inbounds pipeWallMask(fieldNodes[1][i], fieldNodes[3][k]) * rate * (field[i + 1, j, k] - field[i, j, k])
    #on right side of pipe
    else
        return @inbounds pipeWallMask(fieldNodes[1][i], fieldNodes[3][k]) * rate * (field[i - 1, j, k] - field[i, j, k])
    end
end



#TODO: sepparate initial functions and profile functions 



"""NAME OF TRIAL"""
trial_name = "test"



"""SET UP MODEL COMPONENTS"""
#calculating max allowed spacing
surrounding_density_gradient = (getBackgroundDensity(-pipe_top_depth - pipe_length, T_init, S_init) - getBackgroundDensity(-pipe_top_depth, T_init, S_init))/pipe_length #takes average for unaltered water column
surrounding_density_gradient = (getBackgroundDensity(-pipe_top_depth - pipe_length, T_init, S_init_b) - getBackgroundDensity(-pipe_top_depth, T_init, S_init_b))/pipe_length #takes average for unaltered water column
oscillation_angular_frequency =  sqrt((g/1000) * surrounding_density_gradient)
oscillation_period = 2π/oscillation_angular_frequency
@info @sprintf("Buoyancy Oscillation period: %.3f minutes",  oscillation_period/minute)
#grid spacing, NEEDS TO BE LESS THAN 0.01189 for viscosity, but I dont think thats neccessary
max_grid_spacing = 0.05 #manual option
#this sets max grid spacing based off diffusivity through pipe walls, TODO: may be able to increase if we scale thermal diffusivity based off wdith 
#max_grid_spacing = 0.98*sqrt(pipe_data.wall_thermal_diffusivity*0.25*oscillation_period) #TRACER_MIN MAY NEED TO EDIT DEPENDING ON TRACERS: sets max grid spacing to be distance that slowest diffusing tracer travels in 1/4 the oscillation timescale
#max_grid_spacing = 0.98*sqrt(T_diffusivity.molecular*0.25*oscillation_period) #0.0004m RACER_MIN MAY NEED TO EDIT DEPENDING ON TRACERS: sets max grid spacing to be distance that slowest diffusing tracer travels in 1/4 the oscillation timescale
#memory = 64 #computer memory in GB
#max_grid_spacing = max(domain_x, domain_z)/sqrt(((memory/32) * 100000000)) #for 2D, max that computer can handle
# max_grid_spacing = max(domain_x, domain_z)/sqrt(((memory/32) * 100000000)/4) #for 2D, max that computer can handle
#spacing declarations (grid size)
x_grid_spacing = max_grid_spacing;
y_grid_spacing = max_grid_spacing;
z_grid_spacing = max_grid_spacing;
#calculated resolutions (number of grid cells)
x_res = floor(Int, domain_x / x_grid_spacing);
y_res = 4 #4 cells thick in the "flat" direction
z_res = floor(Int, domain_z / z_grid_spacing);
#"flat" dimension domain size
domain_y = y_res * y_grid_spacing
#miscellaneous
pipe_wall_thickness = roundUp(pipe_wall_thickness_intended, x_grid_spacing)
@info @sprintf("Pipe walls are %1.2f meters thick", pipe_wall_thickness)
@info @sprintf("X spacings: %.3e meters | Y spacings: %.3e meters | Z spacings: %.3e meters", x_grid_spacing, y_grid_spacing, z_grid_spacing)
@info @sprintf("X resolution: %.3e | Y resolution: %.3e | Z resolution: %.3e ", x_res, y_res, z_res)
#set up domain grid 
#BOUNDED OPTIONS (LAB_SET)
#3d option (4 cell) 
#domain_grid = RectilinearGrid(CPU(), Float64; size=(x_res, y_res, z_res), x=(0, domain_x), y = (0, domain_y), z=(-domain_z, 0), topology=(Bounded, Periodic, Bounded))
#2d option  
domain_grid = RectilinearGrid(CPU(), Float64; size=(x_res, z_res), x=(0, domain_x), z=(-domain_z, 0), topology=(Bounded, Flat, Bounded))
#PERIODIC OPTION
#domain_grid = RectilinearGrid(CPU(), Float64; size=(x_res, y_res, z_res), x=(0, domain_x), y = (0, domain_y), z=(-domain_z, 0), topology=(Periodic, Periodic, Bounded))
#2d option  
#domain_grid = RectilinearGrid(CPU(), Float64; size=(x_res, z_res), x=(0, domain_x), z=(-domain_z, 0), topology=(Periodic, Flat, Bounded))


#MISCELLANEOUS
clock = Clock{eltype(domain_grid)}(time=0);
timestepper = :QuasiAdamsBashforth2; #default, 3rd order option available 


#ADVECTION SCHEME OPTIONS
#advection = CenteredSecondOrder(); # will be relatively more oscillatory, "overshoots"
#advection = CenteredFourthOrder(); # still overshoots, but should be smoother 
advection = WENO(); #will be smooth, and mimic more diffusion, default is 5th order
#advection = WENO(order = 9)


#BUOYANCY MODEL
# buoyancy = nothing # testing option
 eos = TEOS10EquationOfState()
 buoyancy = SeawaterBuoyancy(equation_of_state=eos)
#buoyancy = SeawaterBuoyancy(equation_of_state = LinearEquationOfState(thermal_expansion = 2e-4, haline_contraction = 78e-5)); #TODO: potentially add more accuracy here, currently set to global average


# #SECTION SETS NEW DIFFUSION VALUES BASED ON A DESIRED GRID SIZE, FOR TESTING PURPOSES, DO NOT USE IF EDDY DIFFUSIVITY IN USE, currently not in use due to no oscillation
# target_spacing = max_grid_spacing
# minimum_diffusivity = (max_grid_spacing^2)/(0.25*oscillation_period) #TODO: check what this is 
# minimum_diffusivity = 1.01*minimum_diffusivity #safety factor 
# #now scale all according to minimum diffusivity 
# #currently setting temperature as standard, and only editing tracers not viscosity, since only temperature goes through tube
# #could also set the smallest diffusivity as the standard, in the case, salinity, may need to edit TRACER_MIN
# reference_diffusivity = T_diffusivity.molecular
# function rescaleDiffusivities(diff, scaleWith, scaleTo)
#     return (diff/scaleWith)*scaleTo
# end
# T_diffusivity.molecular = rescaleDiffusivities(T_diffusivity.molecular, reference_diffusivity, minimum_diffusivity)
# S_diffusivity.molecular = rescaleDiffusivities(S_diffusivity.molecular, reference_diffusivity, minimum_diffusivity)
# pipe_data.wall_thermal_diffusivity = rescaleDiffusivities(pipe_data.wall_thermal_diffusivity, reference_diffusivity, minimum_diffusivity)


#TRACERS & DIFFUSION CLOSURES
tracers = (:T, :S); #temperature, salinity
"""IMPORTANT: order of definition of kappas matter,must define in same order as tracer statement even though finding it later works""";
diffusivity_data = (seawater = seawater_diffusion_data, pipe = pipe_data)
#TODO: this function uses nonlocal data  #working on the named tuple it should take 
#this function also uses nonlocal data for pipe wall thickness and pipe wall thickness intended
function tempDiffusivities(x, y, z, parameters::NamedTuple, wall_indicator::String)
    if (isPipeWall(x, z) || wall_indicator == "WALL")
        return parameters[:pipe].wall_thermal_diffusivity * (pipe_wall_thickness/pipe_wall_thickness_intended) #edit to account for wall thickness
    else 
        return parameters[:seawater].T.molecular
    end
end
tempDiffusivities(x, z, t) = tempDiffusivities(x, 0, z, diffusivity_data, "")
function saltDiffusivities(x, y, z, parameters::NamedTuple)
    if (isPipeWall(x, z))
        return 0
    else 
        return parameters[:seawater].S.molecular
    end
end
saltDiffusivities(x ,z ,t ) = saltDiffusivities(x, 0, z, diffusivity_data)
horizontal_closure = HorizontalScalarDiffusivity(ν=my_viscosity.molecular + my_viscosity.isopycnal, κ=(T=T_diffusivity.molecular + T_diffusivity.isopycnal, S=S_diffusivity.molecular + S_diffusivity.isopycnal)) 
vertical_closure = VerticalScalarDiffusivity(ν=my_viscosity.molecular + my_viscosity.diapycnal, κ=(T=T_diffusivity.molecular + T_diffusivity.diapycnal, S=S_diffusivity.molecular + S_diffusivity.diapycnal)) 
#option for no walls, just molecular diffusivities - TESTED
#closure = ScalarDiffusivity(ν=my_viscosity.molecular, κ=(T=T_diffusivity.molecular, S=S_diffusivity.molecular)) #this made the diffusion time scale way too long
#option for eddy diffusivities 
# closure = (horizontal_closure, vertical_closure)
#option for walls - thermal only
# closure = ScalarDiffusivity(ν=my_viscosity.molecular, κ=(T=tempDiffusivities, S=S_diffusivity.molecular))
#option for walls: thermal and salt
closure = ScalarDiffusivity(ν=my_viscosity.molecular, κ=(T=tempDiffusivities, S=saltDiffusivities))


#BOUNDARY CONDITIONS
#initial gradient dζ/dz, assuming z decreases with depth
#TODO: use this instead of relaxation to imitate infinite reservoir
initial_T_top_gradient = (T_init(0,0) - T_init(0, 0 - z_grid_spacing))/z_grid_spacing
initial_T_bottom_gradient = (T_init(0, domain_z + z_grid_spacing) - T_init(0, domain_z))/z_grid_spacing
initial_S_top_gradient = (S_init(0, 0 - z_grid_spacing) - S_init(0,0))/z_grid_spacing
initial_S_bottom_gradient = (S_init(0, domain_z) - S_init(0, domain_z + z_grid_spacing))/z_grid_spacing
initial_S_top_gradient = (S_init_b(0, 0 - z_grid_spacing) - S_init_b(0,0))/z_grid_spacing
initial_S_bottom_gradient = (S_init_b(0, domain_z) - S_init_b(0, domain_z + z_grid_spacing))/z_grid_spacing

#these two only incoporate constant gradient
# T_bcs = FieldBoundaryConditions(top = GradientBoundaryCondition(initial_T_top_gradient), bottom = GradientBoundaryCondition(initial_T_bottom_gradient))
# S_bcs = FieldBoundaryConditions(top = GradientBoundaryCondition(initial_S_top_gradient), bottom = GradientBoundaryCondition(initial_S_bottom_gradient))
#below incorporates side boundaries as a constant gradient, CANT DO WITH PERIODIC (LAB_SET)
T_bcs = FieldBoundaryConditions(top = GradientBoundaryCondition(initial_T_top_gradient), bottom = GradientBoundaryCondition(initial_T_bottom_gradient), east = ValueBoundaryCondition(T_init_bc), west = ValueBoundaryCondition(T_init_bc))
S_bcs = FieldBoundaryConditions(top = GradientBoundaryCondition(initial_S_top_gradient), bottom = GradientBoundaryCondition(initial_S_bottom_gradient), east = ValueBoundaryCondition(S_init_bc), west = ValueBoundaryCondition(S_init_bc))
boundary_conditions = (T = T_bcs, S = S_bcs)


#MANUALLY SET MAX TIME STEP, FOR FORCING FUNCTIONS RELAXATION RATE 
#this sets it automatically based off predicted minimum delta T, which is predicte via height displacemnt and buoyancy frequency
#TODO: is there a way to making foring timescale a function of model, such that it takes in model's timestep and sets it to be just over that --> look into source code 
v_max_predicted = oscillation_angular_frequency * height_displaced
min_time_step_predicted = (min(x_grid_spacing, z_grid_spacing)*CFL)/v_max_predicted
max_time_step = 2 * min_time_step_predicted#sets max time step to be 2x the predicted minimum time step
#this sets a max ∇T time step manaully.  set this to be larger than what you think the time step will be
#provides the option to use this in determining relaxation rate instead of oscillation time scale --> allows bigger steps than 1/4 oscillation time scale 
# max_time_step = 0.3 #in seconds


#FORCING FUNCTIONS
#option for no forcing
noforcing(x, z, t) = 0
noforcing(x, y, z, t) = noforcing(x, z, t)
#sets velocities inside wall to be 0
max_damping_rate = 1.1*max_time_step #this is the max frequency, eg. can only choose damping timescales bigger  
u_damping_rate = 1/max_damping_rate #relaxes fields on 0.1 second time scale, should be very high
w_damping_rate = 1/max_damping_rate #relaxes fields on 0.11 second time 
u_pipe_wall = Relaxation(rate = u_damping_rate, mask = pipeWallMask)
w_pipe_wall = Relaxation(rate = w_damping_rate, mask = pipeWallMask)
#sets sponge layer for parameters on all water borders
border_damping_rate = 1/max_damping_rate
uw_border = Relaxation(rate = border_damping_rate, mask = waterBorderMask)
T_border = Relaxation(rate = border_damping_rate, mask = waterBorderMask, target = T_init_target)
S_border = Relaxation(rate = border_damping_rate, mask = waterBorderMask, target = S_init_target)
#sets forcing for salinity inside pipe --> NO LONGER USED, USE WALL FORCING INSTEAD, SEE TRACER FORCING IN WALLS BELOW
# pipe_damping_rate = 1/max_damping_rate
# S_pipe = Relaxation(rate = pipe_damping_rate, mask = pipeMask, target = S_init(0, 0, -pipe_bottom_depth))
#sets initial forcing for velocity while pump driven 
pump_info = (Δz = height_displaced, wₚ = initial_pipe_velocity, f_rate = 1/max_damping_rate)
w_pump_forcing = Forcing(w_pump, parameters = pump_info, field_dependencies = :w)
#sets tracer in walls to be forced to that immediately inside pipe
tracer_field_nodes = nodes(domain_grid,(Center(), Center(), Center()); reshape = true)
S_wall_forcing_func(i, j, k, grid, clock, model_fields, rate) = inpenetrable_wall_forcer(x_center, tracer_field_nodes, i, j, k, model_fields.S, rate)
S_wall_forcing = Forcing(S_wall_forcing_func, discrete_form = true, parameters = 1/max_damping_rate)
#forces the outside of the pipe immediately to the surrounding to prevent diffusion
pipe_exterior_damping_rate = 1/max_damping_rate
S_wall_exterior = Relaxation(rate = pipe_exterior_damping_rate, mask = tracerRelaxationMask, target = S_init_target)
#no forcing
# forcing = (u = noforcing,  w = noforcing, T = noforcing, S = noforcing)
# pipe wall velocities only 
# forcing = (u = u_pipe_wall,  w = w_pipe_wall, T = noforcing, S = noforcing)
#pipe wall velocities and property sponge layer 
#forcing = (u = (u_pipe_wall, uw_border),  w = (w_pipe_wall, uw_border), T = T_border, S=S_border)
forcing = (u = (u_pipe_wall, uw_border),  w = (w_pipe_wall, uw_border), T = T_border, S=S_border)
#pipe wall velocities and property sponge layer, and pipe wall relaxation - LATEST TESTED 
forcing = (u = (u_pipe_wall, uw_border),  w = (w_pipe_wall, uw_border), T = T_border, S=(S_border, S_wall_forcing))
#forcing = (u = (u_pipe_wall, uw_border),  w = (w_pipe_wall, uw_border), T = T_border, S=(S_border, S_wall_forcing))
#pipe wall velocities and property sponge layer, pipe wall relaxation, exterior pipe relaxation - CURRENTLY TESTING
#forcing = (u = (u_pipe_wall, uw_border),  w = (w_pipe_wall, uw_border), T = T_border, S=(S_border, S_wall_forcing, S_wall_exterior))
#pipe wall velocities and property sponge layer, and pipe wall relaxation, and pumping
#forcing = (u = (u_pipe_wall, uw_border),  w = (w_pipe_wall, uw_border, w_pump_forcing), T = T_border, S=(S_border, S_wall_forcing))
#pipe wall velocities and property sponge layer, and internal pipe relaxation
#forcing = (u = (u_pipe_wall, uw_border),  w = (w_pipe_wall, uw_border), T = T_border, S=(S_border, S_pipe))
#pipe wall velocities and property sponge layer, and internal pipe relaxation, and pumping
#forcing = (u = (u_pipe_wall, uw_border),  w = (w_pipe_wall, uw_border, w_pump_forcing), T = T_border, S=(S_border, S_pipe))
#pipe wall velocities, property sponge layer, internal pipe relaxation & initial pump velocity 
#forcing = (u = (u_pipe_wall, uw_border),  w = (w_pipe_wall, uw_border, w_pump_forcing), T = T_border, S=(S_border, S_pipe))
#BIOGEOCHEMISTRY
# #biogeochemistry =  LOBSTER(; domain_grid); #not yet used at all



#sets up model
model = NonhydrostaticModel(; grid=domain_grid, clock, advection, buoyancy, tracers, timestepper, closure, forcing, boundary_conditions)
@info "model made"
#helper functions
density_operation = seawater_density(model; geopotential_height)


set!(model, T=T_init, S=S_init, w=w_init)
ρ_initial = Field(density_operation)
compute!(ρ_initial)
@info "initial conditions set"


#setting time steps 
min_grid_spacing = min(minimum_xspacing(model.grid), minimum_zspacing(model.grid))
initial_travel_velocity = initial_pipe_velocity # a number just to set an initial time step, determined from initial runs starting around 1ms, was 0.01 first 
initial_advection_time_scale = min_grid_spacing/initial_travel_velocity
# diffusion_time_scale = (min_grid_spacing^2)/model.closure.κ.T #TRACER_MIN, set to tracer with biggest kappa
# diffusion_time_scale = (min_grid_spacing^2)/model.closure[1].κ.T #TRACER_MIN, set to tracer with biggest kappa, to use when horizontal and vertical diffusivities are used
diffusion_time_scale = (min_grid_spacing^2)/model.closure.κ.T(0, 0, 0, diffusivity_data, "WALL") #TRACER_MIN, set to tracer with biggest kappa, to use when using function based diffusivity
viscous_time_scale = (min_grid_spacing^2)/model.closure.ν
# surrounding_density_gradient = (getBackgroundDensity(-pipe_top_depth - pipe_length, T_init, S_init) - getBackgroundDensity(-pipe_top_depth, T_init, S_init))/pipe_length#takes average for unaltered water column
# initial_pipe_density = getMaskedAverageTracer(pipeMask, ρ_initial)
# initial_oscillation_period = 2π/sqrt((g/initial_pipe_density) * surrounding_density_gradient) #this is more accurate than it needs to be, can replace initial pipe density with 1000


#for intial time step to be ok, next step need to increase, currently its 500ms, seems to increase)
initial_time_step = 0.5*min(0.2 * min(diffusion_time_scale, oscillation_period, viscous_time_scale, initial_advection_time_scale), max_time_step)
"""IMPORTANT, diffusion cfl does not work with functional kappas, need to manually set max step"""
new_max_time_step = min(0.2 * diffusion_time_scale, max_time_step) #TRACER_MIN, uses a cfl of 0.2, put in diffusion time scale of tracer with biggest kappa, or viscosity
simulation_duration = 1day
run_duration = 6hour
run_duration = 10minute


#running model
simulation = Simulation(model, Δt=initial_time_step, stop_time=simulation_duration, wall_time_limit=run_duration) # make initial delta t bigger
timeWizard = TimeStepWizard(cfl=CFL, diffusive_cfl = CFL, max_Δt = new_max_time_step) #TODO: set max delta t?
simulation.callbacks[:timeWizard] = Callback(timeWizard, IterationInterval(4))
progress_message(sim) = @printf("Iteration: %04d, time: %s, Δt: %s, wall time: %s\n", iteration(sim), prettytime(sim), prettytime(sim.Δt), prettytime(sim.run_wall_time))
add_callback!(simulation, progress_message, IterationInterval(50))

#fields = Dict("u" => model.velocities.u, "w" => model.velocities.w, "T" => model.tracers.T, "S" => model.tracers.S)
w = model.velocities.w;
u = model.velocities.u;
T = model.tracers.T;
S = model.tracers.S;
ζ = Field(-∂x(w) + ∂z(u)) #vorticity in y 
ρ = Field(density_operation)
filename = joinpath("Trials",trial_name)
#simulation.output_writers[:outputs] = JLD2OutputWriter(model, (; u, w, T, S, ζ, ρ); filename, schedule=IterationInterval(10), overwrite_existing=true) 
#time interval option
simulation.output_writers[:outputs] = JLD2OutputWriter(model, (; u, w, T, S, ζ, ρ); filename, schedule=TimeInterval(10), overwrite_existing=true) #can also set to TimeInterval
#time average perhaps?
run!(simulation; pickup=false)



#visualize simulation
output_filename = filename * ".jld2"
u_t = FieldTimeSeries(output_filename, "u")
w_t = FieldTimeSeries(output_filename, "w")
ζ_t = FieldTimeSeries(output_filename, "ζ")
T_t = FieldTimeSeries(output_filename, "T")
S_t = FieldTimeSeries(output_filename, "S")
ρ_t = FieldTimeSeries(output_filename, "ρ")
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
@info "finished extracting data as arrays"
#how much of data set to plot 
plot_until_time = -1 #enter -1 to plot whole thing
lastFrame = getIndexFromTime(plot_until_time, times)
num_Data_Points = lastFrame
#very inefficient way of getting max/min, need to update
T_range = getMaxAndMin(num_Data_Points, T_t)
S_range = getMaxAndMin(num_Data_Points, S_t)
u_range = getMaxAndMin(num_Data_Points, u_t)
w_range = getMaxAndMin(num_Data_Points, w_t)
ρ_range = getMaxAndMin(num_Data_Points, ρ_t)
ζ_range = getMaxAndMin(num_Data_Points, ζ_t)
@info "finished getting max and min of each"


#making animations 
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
frames = 1:lastFrame
record(fig, filename * "properties.mp4", frames, framerate=8) do i
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
frames = 1:lastFrame
record(fig, filename * "velocities.mp4", frames, framerate=8) do i
    @info string("Plotting frame ", i, " of ", frames[end])
    n[] = i
end
@info @sprintf("Done. Simulation total run time %.3f seconds | %.3f minutes", times[end], times[end]/minute)

























# xT, yT, zT = nodes(T_t[1])
# ax_T = Axis(fig[3, 1]; title = "temperature", axis_kwargs...)
# hm_T = heatmap!(ax_T, xT, zT, Tₙ; colorrange = T_range, colormap = :thermal)
# Colorbar(fig[3,2], hm_T, label = "C")
#TODO: visualize streamlines & how they change
#visualize distributions
density_field = Field(density_operation)
compute!(density_field)
#function takes a model
function plotInitialConditions(m)
    geopotential_height = 0
    density_operation = density_operation = seawater_density(m; geopotential_height)
    density_field = Field(density_operation)
    compute!(density_field)
    fig = Figure()
    title = "initial conditions"
    Label(fig[0, :], title)
    x = xnodes(m.tracers.T)
    z = znodes(m.tracers.T)
    axis_kwargs = (xlabel="x (m)", ylabel="z (m)", width=400)
    ax1 = Axis(fig[1, 1]; title="Temperature (C)", axis_kwargs...)
    T_init_matrix = interior(m.tracers.T, :, 1, :)
    T_lims = (minimum(T_init_matrix), maximum(T_init_matrix));
    hm1 = heatmap!(ax1, x, z, T_init_matrix, colorrange=T_lims, colormap=:thermal, interpolate=true)
    Colorbar(fig[1, 2], hm1)
    ax2 = Axis(fig[2, 1]; title="Salinity (ppt)", axis_kwargs...)
    S_init_matrix = interior(m.tracers.S, :, 1, :)
    S_lims = (minimum(S_init_matrix), maximum(S_init_matrix));
    hm2 = heatmap!(ax2, x, z, S_init_matrix, colorrange=S_lims, colormap=:haline, interpolate=true)
    Colorbar(fig[2, 2], hm2)
    ax3 = Axis(fig[3, 1]; title="Density (kg/m^3)", axis_kwargs...)
    ρ_init_matrix = interior(density_field, :, 1, :)
    ρ_lims = (minimum(ρ_init_matrix), maximum(ρ_init_matrix));
    hm3 = heatmap!(ax3, x, z, ρ_init_matrix, colorrange=ρ_lims, colormap=Reverse(:viridis), interpolate=true)
    Colorbar(fig[3, 2], hm3)
    return fig
end 
function plotTracerMask(mask::Function, m)
    fig = Figure()
    title = "tracer mask"
    Label(fig[0, :], title)
    x = xnodes(m.tracers.T)
    z = znodes(m.tracers.T)
    myNodes = nodes(model.grid, (Center(), Center(), Center()), reshape = true)
    maskArray = zeros(Float64, size(myNodes[1])[1], size(myNodes[3])[3])
    for i in eachindex(myNodes[1])
        for j in eachindex(myNodes[3])
            maskArray[i, j] = mask(myNodes[1][i], myNodes[3][j])
        end
    end
    axis_kwargs = (xlabel="x (m)", ylabel="z (m)", width=400)
    ax = Axis(fig[1, 1]; title="Temperature (C)", axis_kwargs...)
    lims = (minimum(maskArray), maximum(maskArray));
    hm = heatmap!(ax, x, z, maskArray, colorrange=lims, colormap=:grays, interpolate=true)
    Colorbar(fig[1, 2], hm)
    return fig
end
#check that model has reset - plot velocities
fig = Figure()
title = "initial velocities"
Label(fig[0, :], title)
x = xnodes(model.tracers.T)
z = znodes(model.tracers.T)
ax1 = Axis(fig[2, 1]; title="u velocity (m/s)", axis_kwargs...)
u_init_matrix = interior(model.velocities.u, :, 1, :)
u_lims = (minimum(u_init_matrix), maximum(u_init_matrix));
hm1 = heatmap!(ax1, x, z, u_init_matrix, colorrange=u_lims, colormap=Reverse(:balance), interpolate=true)
Colorbar(fig[2, 2], hm1)
axis_kwargs = (xlabel="x (m)", ylabel="z (m)", width=400)
ax2 = Axis(fig[1, 1]; title="w velocity (m/s)", axis_kwargs...)
w_init_matrix = interior(model.velocities.w, :, 1, :)
w_lims = (minimum(w_init_matrix), maximum(w_init_matrix));
hm2 = heatmap!(ax2, x, z, w_init_matrix, colorrange=w_lims, colormap=Reverse(:balance), interpolate=true)
Colorbar(fig[1, 2], hm2)
display(fig)
#plotting streamlines: need to write up an troubleshoot 
x_spacings = xspacings(domain_grid, Center(), Center(), Center())
z_spacings = zspacings(domain_grid, Center(), Center(), Center())
u_init_matrix = interior(model.velocities.u, :, 1, :)
w_init_matrix = interior(model.velocities.w, :, 1, :)
x = xnodes(model.tracers.T)
z = znodes(model.tracers.T)
function streamFunc(xCoord, zCoord)
    #need to add interpolation 
    x_ind = floor(Int, (floor(Int, 10xCoord) / 10 + x_spacings) / x_spacings)
    z_ind = -floor(Int, (floor(Int, 10zCoord) / 10 + z_spacings) / z_spacings)
    newX = xCoord + 0.01 * (u_init_matrix[x_ind, z_ind])
    newZ = zCoord + 0.01 * (w_init_matrix[x_ind, z_ind])
    return Point2f(newX, newZ)
end
fig = Figure(size=(300, 300))
title = "streamlines"
Label(fig[0, :], title)
ax = Axis(fig[1, 1]; title="streamlines", axis_kwargs...)
streamplot(streamFunc, 0 .. domain_x, -domain_z .. 0, colormap=:magma)
fig
#TODO: why does it seem to get so high, plot total heat to see if conservative 
#testing
#setting up pipe
# domain_arr = nodes(domain_grid, locs; reshape = true)
# pipe_interior_grid_x = AbstractArray{Float64, 3}
x_spacings = xspacings(domain_grid, Center(), Center(), Center())
z_spacings = zspacings(domain_grid, Center(), Center(), Center())