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

#IMPORTANT, IF WRITE DISCRETE FORCER HERE, NEED TO BE CAREFUL ABOUT CUARRAY VS ARRAY AND ADAPT THE CUARRAY OVER
"""NAME"""
trial_name = "gpu zero side viscosity double cell model"

"""CONSTANTS AND PHYSICAL PROPERTIES"""
#utility constants
const day = 86400;
const hour = 3600;
const minute = 60;
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
const sw_viscosity_molecular = 1.05e-6
const sw_T_diffusivity_molecular = 1.46e-7
const sw_S_diffusivity_molecular = 1.3e-9
const sw_diffusivity_data = SeawaterDiffusivityData(sw_viscosity_molecular, sw_T_diffusivity_molecular, sw_S_diffusivity_molecular)

"""DOMAIN SIZE & SETUP"""
const domain_x = 20;
const domain_z = 20;
const x_center = domain_x/2;

"""PIPE SIZE AND SETUP"""
struct PipeWallData
    thermal_diffusivity :: Float64
    thickness :: Float64
end
const pipe_radius = 0.25
const pipe_length = 10
const pipe_top_depth = 4
const pipe_wall_thickness_intended = 0.01
const pipe_bottom_depth = pipe_top_depth + pipe_length
const wall_material_ρ = 8900
const wall_material_cₚ = 376.812 
const wall_material_k = 50.208
const pipe_data = PipeWallData(wall_material_k/(wall_material_cₚ*wall_material_ρ), pipe_wall_thickness_intended)

"""SET BACKGROUND TEMPERATURE AND SALINITY GRADIENTS"""
const T_top = 21.67;
const T_bot = 11.86;
const S_bot = 34.18;
const S_top = 35.22;
const delta_z = 200; 

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
        return S_top - ((S_bot - S_top) / delta_z)*(-pipe_top_depth)
    elseif (z < -pipe_bottom_depth)
        return S_top - ((S_bot - S_top) / delta_z)*(-pipe_bottom_depth)
    else
        return S_top - ((S_bot - S_top) / delta_z)z
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
function getMaxAndMin(numPoints, dataSeries)
    myMax = maximum(interior(dataSeries[1], :, 1, :))
    myMin = minimum(interior(dataSeries[1], :, 1, :))
    for i in 2:numPoints
        myMax = max(maximum(interior(dataSeries[i], :, 1, :)),myMax )
        myMin = min(minimum(interior(dataSeries[i], :, 1, :)), myMin)
    end
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
function getBackgroundDensity(z, Tfunc::Function, Sfunc::Function)#takes a positive depth
    eos = TEOS10EquationOfState() 
    S_abs = gsw_sa_from_sp(Sfunc(0,z), gsw_p_from_z(z, 30), 31, -30) #random lat and long in ocean
    T_consv = gsw_ct_from_t(S_abs, Tfunc(0,z), gsw_p_from_z(z, 30)) #set to 30 degrees north
    return TEOS10.ρ(T_consv, S_abs, z, eos) #note that this uses insitu s and T instead of conservative and absolute, which is what the function calls for 
end
function checkMemory(capacity, x_cells, z_cells)
    numCells = x_cells * z_cells
    possibleCells = (capacity/32) * 1000000 
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
@info @sprintf("Buoyancy Oscillation period: %.3f minutes",  oscillation_period/minute)
#grid spacing
const max_grid_spacing = 0.05 #TODO: figure out what this needs to be set to 
const x_grid_spacing = max_grid_spacing;
const z_grid_spacing = 2*max_grid_spacing;
#resolution
x_res = floor(Int, domain_x / x_grid_spacing);
z_res = floor(Int, domain_z / z_grid_spacing);
#call error if it is too large for model
const GPU_memory = 12 #GB
checkMemory(GPU_memory, x_res, z_res)
#pipe wall thickness, as set by the model
const pipe_wall_thickness = roundUp(pipe_wall_thickness_intended, x_grid_spacing)
@info @sprintf("Pipe walls are %1.2f meters thick", pipe_wall_thickness)
#time stepping 


#location with equivalent density
#TODO: potentially fix that height displaced with an estimation of where the density is equal 
# v_max_predicted = oscillation_angular_frequency * height_displaced
v_max_predicted = 0.025
min_time_step_predicted = (min(x_grid_spacing, z_grid_spacing)*CFL)/v_max_predicted
max_time_step = 2 * min_time_step_predicted
#damping rate for forcers  
min_relaxation_timescale = 1.1 * max_time_step
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
function pipeWallMask(x,z)
    if (isPipeWall(x, z))
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
tracerRelaxationMaskDomainTwo(x, z) = tracerRelaxationMaskDomainTwo(x, 0, z)
velocityRelaxationMaskDomainOne(x, z) = velocityRelaxationMaskDomainOne(x, 0, z)
tracerRelaxationMaskDomainThree(x, z) = tracerRelaxationMaskDomainThree(x, 0, z)

"""INITIAL CONDITIONS"""
function T_init(x, y, z)
    return TWaterColumn(z)
end

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
T_init(x, z) = T_init(x, 0, z)
S_init(x, z) = S_init(x, 0, z)
A_init(x, z) = A_init(x, 0, z)

"""MODEL BOUNDARY CONDITIONS SETUP"""
#get gradients for boundary conditons
initial_T_top_gradient = (TWaterColumn(0) - TWaterColumn(0 - z_grid_spacing))/z_grid_spacing
initial_T_bottom_gradient = (TWaterColumn(domain_z + z_grid_spacing) - TWaterColumn(0domain_z))/z_grid_spacing
initial_S_top_gradient = (SWaterColumn(0, 0 - z_grid_spacing) - SWaterColumn(0,0))/z_grid_spacing
initial_S_bottom_gradient = (SWaterColumn(0, domain_z) - SWaterColumn(0, domain_z + z_grid_spacing))/z_grid_spacing

"""FORCING FUNCTIONS"""
#none currently used 
TWaterColumnTarget(x, z, t) = TWaterColumn(z)
SWaterColumnTarget(x, z, t) = SWaterColumn(z)

"""MODEL DIFFUSIVITY SETUP"""
const diffusivity_data = (seawater = sw_diffusivity_data, pipe = pipe_data)
const pipeWallThickness = (actual = pipe_wall_thickness, intended = pipe_wall_thickness_intended)
function tempDiffusivities(x, y, z, diffusivities::NamedTuple, wallThickness::NamedTuple, wall_indicator::String)
    if (isPipeWall(x, z) || wall_indicator == "WALL")
        return diffusivities[:pipe].thermal_diffusivity * (wallThickness[:actual]/wallThickness[:intended]) #edit to account for wall thickness
    else 
        return diffusivities[:seawater].T
    end
end
function saltDiffusivities(x, y, z, diffusivities::NamedTuple)
    if (isPipeWall(x, z))
        return 0
    else 
        return diffusivities[:seawater].S
    end
end
function myViscosity(x, y, z, diffusivities::NamedTuple)
    if (isPipeWall(x, z))
        return 0
    elseif (isSidePipe(x, z))
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
buoyancy = SeawaterBuoyancy(equation_of_state=eos)

tracers = (:T, :S, :A)
closure = ScalarDiffusivity(ν=myViscosity, κ=(T=tempDiffusivities, S=saltDiffusivities, A = 0))

T_bcs = FieldBoundaryConditions(top = GradientBoundaryCondition(initial_T_top_gradient), bottom = GradientBoundaryCondition(initial_T_bottom_gradient))
S_bcs = FieldBoundaryConditions(top = GradientBoundaryCondition(initial_S_top_gradient), bottom = GradientBoundaryCondition(initial_S_bottom_gradient))
boundary_conditions = (T = T_bcs, S = S_bcs)

#TODO: combine these when set 
uw_pipe_wall_forcer = Relaxation(rate = max_relaxation_rate, mask = pipeWallMask)
T_domain_forcer = Relaxation(rate = max_relaxation_rate, mask = tracerRelaxationMaskDomainTwo, target = TWaterColumnTarget)
S_domain_forcer = Relaxation(rate = max_relaxation_rate, mask = tracerRelaxationMaskDomainTwo, target = SWaterColumnTarget)
S_side_forcer = Relaxation(rate = max_relaxation_rate, mask = sidePipeMask, target = SWaterColumn(-pipe_top_depth))
uw_domain_forcer = Relaxation(rate = max_relaxation_rate, mask = velocityRelaxationMaskDomainOne)
forcing = (u = uw_pipe_wall_forcer,  w = uw_pipe_wall_forcer, T = (T_domain_forcer), S = (S_domain_forcer, S_side_forcer))

model = NonhydrostaticModel(; grid=domain_grid, clock, advection, buoyancy, tracers, timestepper, closure, forcing, boundary_conditions)
density_operation = seawater_density(model; geopotential_height)
@info "model made"

"""SET UP INITIAL CONDITIONS"""
set!(model, T= T_init, S=S_init, A = A_init)
set!(model, T= T_init, S=S_init, A = A_init) #ASK  - why set twice ??
@info "initial conditions set"


"""SETTING UP SIMULATION"""
#how long you want to run simulation for
simulation_duration = 1day
run_duration = 2hour
output_interval = 1minute

#finding various time scales 
min_grid_spacing = min(minimum_xspacing(model.grid), minimum_zspacing(model.grid))
initial_travel_velocity_fake = 0.02 # a number just to set an initial time step, i set it to be around max for safety
initial_advection_time_scale = min_grid_spacing/initial_travel_velocity_fake
diffusion_time_scale = (min_grid_spacing^2)/model.closure.κ.T(0, 0, 0, diffusivity_data, pipeWallThickness, "WALL") #TRACER_MIN, set to tracer with biggest kappa, to use when using function based diffusivity
viscous_time_scale = (min_grid_spacing^2)/model.closure.ν(0, 0, 0)
initial_time_step = 0.5*min(0.2 * min(diffusion_time_scale, oscillation_period, viscous_time_scale, initial_advection_time_scale), max_time_step)
"""IMPORTANT, diffusion cfl does not work with functional kappas, need to manually set max step"""
new_max_time_step = min(0.2 * diffusion_time_scale, max_time_step) #TRACER_MIN, uses a cfl of 0.2, put in diffusion time scale of tracer with biggest kappa, or viscosity

#set up simulation & timewizard
simulation = Simulation(model, Δt=initial_time_step, stop_time=simulation_duration, wall_time_limit=run_duration) # make initial delta t bigger

#various callbacks
#timewizard
timeWizard = TimeStepWizard(cfl=CFL, diffusive_cfl = CFL, max_Δt = new_max_time_step, max_change = 1.1, min_change = 0.5) 
simulation.callbacks[:timeWizard] = Callback(timeWizard, IterationInterval(4))
#progress Message
progress_message(sim) = @printf("Iteration: %04d, time: %s, Δt: %s, wall time: %s\n", iteration(sim), prettytime(sim), prettytime(sim.Δt), prettytime(sim.run_wall_time))
add_callback!(simulation, progress_message, IterationInterval(50))
#set up output_writers
w = model.velocities.w;
u = model.velocities.u;
T = model.tracers.T;
S = model.tracers.S;
ζ = Field(-∂x(w) + ∂z(u)) #vorticity in y 
ρ = Field(density_operation)
A = model.tracers.A
filename = joinpath("Trials",trial_name)
#simulation.output_writers[:outputs] = JLD2OutputWriter(model, (; u, w, T, S, ζ, ρ); filename, schedule=IterationInterval(10), overwrite_existing=true) 
#time interval option
simulation.output_writers[:outputs] = JLD2OutputWriter(model, (; u, w, T, S, ζ, ρ, A); filename, schedule=TimeInterval(output_interval), overwrite_existing=true) #can also set to TimeInterval

#run simulation
run!(simulation; pickup=false)


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
@info "finished extracting data as arrays"
#how much of data set to plot 
num_Data_Points = length(times)
#very inefficient way of getting max/min, need to update
T_range = getMaxAndMin(num_Data_Points, T_t)
S_range = getMaxAndMin(num_Data_Points, S_t)
u_range = getMaxAndMin(num_Data_Points, u_t)
w_range = getMaxAndMin(num_Data_Points, w_t)
ρ_range = getMaxAndMin(num_Data_Points, ρ_t)
ζ_range = getMaxAndMin(num_Data_Points, ζ_t)
A_range = getMaxAndMin(num_Data_Points, A_t)
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
frames = 1:length(times)
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
frames = 1:length(times)
record(fig, filename * "velocities.mp4", frames, framerate=8) do i
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
record(fig, filename * "tracer.mp4", frames, framerate=8) do i
    @info string("Plotting frame ", i, " of ", frames[end])
    n[] = i
end
@info @sprintf("Done. Simulation total run time %.3f seconds | %.3f minutes", times[end], times[end]/minute)

































"""MORE PLOTTING"""
function plotInitialConditions(m)
    geopotential_height = 0
    density_operation = density_operation = seawater_density(m; geopotential_height)
    density_field = Field(density_operation)
    compute!(density_field)
    fig = Figure(size=(600, 800))
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
    ax4 = Axis(fig[4, 1]; title="Tracer", axis_kwargs...)
    A_init_matrix = interior(m.tracers.A, :, 1, :)
    A_lims = (minimum(A_init_matrix), maximum(A_init_matrix));
    hm4 = heatmap!(ax4, x, z, A_init_matrix, colorrange=A_lims, colormap=:matter, interpolate=true)
    Colorbar(fig[4, 2], hm4)
    return fig
end 
function plotTracerMask(mask::Function, m)
    fig = Figure()
    title = "tracer mask for "*String(Symbol(mask))
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
    ax = Axis(fig[1, 1]; title="Mask Intensity", axis_kwargs...)
    lims = (minimum(maskArray), maximum(maskArray));
    hm = heatmap!(ax, x, z, maskArray, colorrange=lims, colormap=:grays, interpolate=true)
    Colorbar(fig[1, 2], hm)
    return fig
end

plotInitialConditions(model)


#=notes --> another way to call functions, is to have one then call another that 
passes the "global" parameters in  
for example
T = 3
function myFunc(x, VAR)
    return x*VAR
end
myFunc(x) = myFUNC(x, T)
=#
