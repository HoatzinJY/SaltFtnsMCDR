using Pkg
using Oceananigans
using OceanBioME
using TimesDates
using CairoMakie
using Printf
using Oceananigans, SeawaterPolynomials.TEOS10
using Oceananigans.Models: seawater_density

#TODO: organizize file, declare variables as global and make function file 
#constants and operations
#TODO: figure out why it seems to be diverging, timestep too large, or, diverges over approx 50min
#linked with above, dampening?, should besomething to do with viscosity? 
const day = 86400;
const hour = 3600;
const minute = 60;
geopotential_height = 0; # sea surface height for potential density calculations
roundUp(num::Float64, base) = ceil(Int, num / base) * base


#resolution parameters TODO determine best ones
#when set to 0.1, unstable, seems to be increasing velocity drastically
#note: be careful of time resolution needed, spacing resoslution needed
#user set
x_grid_spacing = 0.1;# in meters
z_grid_spacing = 0.1;# in meters
domain_x = 15; # width, in meters
domain_z = 15; #height, in meters
initial_time_step = 1;
max_time_step = 1minute;
simulation_duration = 1hour;
run_duration = 1minute;

#pipe parameters
pipe_radius = 0.5;
pipe_length = 10;
pipe_top_depth = 1;
pipe_wall_thickness_intended = 0.01;
pipe_wall_thickness = roundUp(pipe_wall_thickness_intended, x_grid_spacing)
@info @sprintf("Pipe walls are %1.2f meters thick", pipe_wall_thickness)

#pumping parameters
height_displaced = 3;
initial_pipe_velocity = 0.001; #not yet used 

#calculated
x_res = floor(Int, domain_x / x_grid_spacing);
z_res = floor(Int, domain_z / z_grid_spacing);
locs = (Center(), Face(), Center());
x_center = domain_x / 2;
function isInsidePipe(x, z)
    if (x > (x_center - pipe_radius) && x < (x_center + pipe_radius) && z < -pipe_top_depth && z > -(pipe_top_depth + pipe_length))
        return true
    else
        return false
    end
end
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
function pipeWallMask(x,y,z)
    if (isPipeWall(x, z))
        return 1
    else
        return 0
    end
end


#setting up model components
domain_grid = RectilinearGrid(CPU(), Float64; size=(x_res, z_res), x=(0, domain_x), z=(-domain_z, 0), topology=(Bounded, Flat, Bounded))
clock = Clock{eltype(domain_grid)}(time=0);
advection = CenteredSecondOrder(); #default, not sure which one to choose
eos = TEOS10EquationOfState()
buoyancy = SeawaterBuoyancy(equation_of_state=eos)
#buoyancy = SeawaterBuoyancy(equation_of_state = LinearEquationOfState(thermal_expansion = 2e-4, haline_contraction = 78e-5)); #TODO: potentially add more accuracy here, currently set to global average
tracers = (:T, :S); #temperature, salinity
timestepper = :QuasiAdamsBashforth2; #default, 3rd order option available 
closure = ScalarDiffusivity(ν=1e-6, κ=(S=1e-7, T=1.45e-10));
# #biogeochemistry =  LOBSTER(; domain_grid); #not yet used at all
# background_fields::water_column_profiles = water_column_profiles(); #TODO, add in background t & s profiles to serve as a "contant" beyond pipe 
# forcing = wallForcers(); #for imaginary pipe walls
#TODO: more accurate closure/diffusivity, set constants as functions
#uses relaxation to set pipe wall velocities to 0 TODO: look into implementation of relaxation
u_damping_rate = 1/0.1 #relaxes fields on 0.1 second time scale
w_damping_rate = 1/1 #relaxes fields on 1 second time scale
u_pipe_wall = Relaxation(rate = u_damping_rate, mask = pipeWallMask)
w_pipe_wall = Relaxation(rate = w_damping_rate, mask = pipeWallMask)
forcing = (u = u_pipe_wall, w = w_pipe_wall)

#sets up model
#following are considered negligable/not accounted for: coriolis, stokes drift
#following are determinined automatically: pressure solver 
#have yet to think about following: auxiliary fields

model = NonhydrostaticModel(; grid=domain_grid, clock, advection, buoyancy = nothing, tracers, timestepper, closure, forcing)
@info "model made"
#helper functions
#density_operation = seawater_density(model; geopotential_height)
#rounds number num up to nearest multiple of base
#returns true or false for if coordinate is inside pipe, relies on center grid 

#setting up just a basic gradient for water column 
T_top = 21.67;
T_bot = 11.86;
S_bot = 35.22;
S_top = 34.18;
delta_z = 200;
#sets initial t and s
#TODO:set up container functions for these to allow moving to another file 
function T_init(x, z)
    #inside pipe
    if (isInsidePipe(x, z))
        return T_top - ((T_bot - T_top) / delta_z) * (z + height_displaced)
        #outside pipe
    else
        return T_top - ((T_bot - T_top) / delta_z)z
    end
end
function S_init(x, z)
    if (isInsidePipe(x, z))
        return S_top - ((S_bot - S_top) / delta_z) * (z + height_displaced)
    else
        return S_top - ((S_bot - S_top) / delta_z)z
    end
end
function w_init(x, z)
    if (isInsidePipe(x, z))
        #return initial_pipe_velocity;
        return 0
    else
        return 0
    end
end
set!(model, T=T_init, S=S_init, w=w_init)
@info "initial ocnditions set"

#running model
simulation = Simulation(model, Δt=initial_time_step, stop_time=simulation_duration, wall_time_limit=run_duration) # make initial delta t bigger
timeWizard = TimeStepWizard(cfl=0.33, max_Δt=max_time_step) #TODO: set max delta t?
simulation.callbacks[:timeWizard] = Callback(timeWizard, IterationInterval(4))
progress_message(sim) = @printf("Iteration: %04d, time: %s, Δt: %s, wall time: %s\n", iteration(sim), prettytime(sim), prettytime(sim.Δt), prettytime(sim.run_wall_time))
add_callback!(simulation, progress_message, IterationInterval(10))
#log progress --> TODO: add this

#fields = Dict("u" => model.velocities.u, "w" => model.velocities.w, "T" => model.tracers.T, "S" => model.tracers.S)
w = model.velocities.w;
u = model.velocities.u;
T = model.tracers.T;
S = model.tracers.S;
ζ = Field(-∂x(w) + ∂z(u)) #vorticity in y 
#ρ = Field(density_operation)
filename = "diffusiontest2"
simulation.output_writers[:outputs] = JLD2OutputWriter(model, (; u, w, T, S, ζ); filename, schedule=TimeInterval(10), overwrite_existing=true)
#time average perhaps?

run!(simulation; pickup=false)





#visualize simulation
#may want to average

output_filename = filename * ".jld2"
u_t = FieldTimeSeries(output_filename, "u")
w_t = FieldTimeSeries(output_filename, "w")
ζ_t = FieldTimeSeries(output_filename, "ζ")
T_t = FieldTimeSeries(output_filename, "T")
S_t = FieldTimeSeries(output_filename, "S")
# ρ_t = FieldTimeSeries(output_filename, "ρ")
times = u_t.times
n = Observable(1)
Observable(1)

uₙ = @lift interior(u_t[$n], :, 1, :)
wₙ = @lift interior(w_t[$n], :, 1, :)
ζₙ = @lift interior(ζ_t[$n], :, 1, :)
Tₙ = @lift interior(T_t[$n], :, 1, :)
Sₙ = @lift interior(S_t[$n], :, 1, :)
# ρₙ = @lift interior(ρ_t[$n], :, 1, :)

num_Data_Points = length(times)
#very inefficient way of getting max/min, need to update
function getMaxAndMin(numPoints, dataSeries)
    max = maximum(dataSeries[1])
    min = minimum(dataSeries[1])
    for i in 2:numPoints
        myMax = maximum(dataSeries[i])
        myMin = minimum(dataSeries[i])
        if (myMax > max)
            max = myMax
        end
        if (myMin < min)
            min = myMin
        end
    end
    return (min, max)
end
T_range = getMaxAndMin(num_Data_Points, T_t)
S_range = getMaxAndMin(num_Data_Points, S_t)
u_range = getMaxAndMin(num_Data_Points, u_t)
w_range = getMaxAndMin(num_Data_Points, w_t)
# ρ_range = getMaxAndMin(num_Data_Points, ρ_t)
ζ_range = getMaxAndMin(num_Data_Points, ζ_t)









#making animations 
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
plot_until_time = -1 #enter -1 to plot whole thing
lastFrame = getIndexFromTime(plot_until_time, times)

#velocities
fig = Figure(size=(600, 900))
title = @lift @sprintf("t = %1.2f minutes", round(times[$n] / minute, digits=2))
axis_kwargs = (xlabel="x (m)", ylabel="z (m)", width=400)
fig[1, :] = Label(fig, title)

xw, yw, zw = nodes(w_t[1])
ax_w = Axis(fig[2, 1]; title="w velocity", axis_kwargs...)
hm_w = heatmap!(ax_w, xw, zw, wₙ; colorrange=w_range, colormap=:balance) #note that this is still using old grid from T, S, initial, may need to recompute x and z using specific nodes 
Colorbar(fig[2, 2], hm_w, label="m/s")

xu, yu, zu = nodes(u_t[1])
ax_u = Axis(fig[3, 1]; title="u velocity", axis_kwargs...)
hm_u = heatmap!(ax_u, xu, zu, uₙ; colorrange=u_range, colormap=:balance) #note that this is still using old grid from T, S, initial, may need to recompute x and z using specific nodes 
Colorbar(fig[3, 2], hm_u, label="m/s")

xζ, yζ, zζ = nodes(ζ_t[1])
ax_ζ = Axis(fig[4, 1]; title="vorticity", axis_kwargs...)
hm_ζ = heatmap!(ax_ζ, xζ, zζ, ζₙ; colorrange=ζ_range, colormap=:balance) #note that this is still using old grid from T, S, initial, may need to recompute x and z using specific nodes 
Colorbar(fig[4, 2], hm_ζ, label="rot/sec")

fig

@info "Making velocities animation from data"
frames = 1:lastFrame
record(fig, filename * "velocities.mp4", frames, framerate=8) do i
    @info string("Plotting frame ", i, " of ", frames[end])
    n[] = i
end

#properties
fig = Figure(size=(600, 900))
title = @lift @sprintf("t = %1.2f minutes", round(times[$n] / minute, digits=2))
axis_kwargs = (xlabel="x (m)", ylabel="z (m)", width=400)
fig[1, :] = Label(fig, title)

xT, yT, zT = nodes(T_t[1])
ax_T = Axis(fig[2, 1]; title="temperature", axis_kwargs...)
hm_T = heatmap!(ax_T, xT, zT, Tₙ; colorrange=T_range, colormap=:thermal) #note that this is still using old grid from T, S, initial, may need to recompute x and z using specific nodes 
Colorbar(fig[2, 2], hm_T, label="C")

xS, yS, zS = nodes(S_t[1])
ax_S = Axis(fig[3, 1]; title="salinity", axis_kwargs...)
hm_S = heatmap!(ax_S, xS, zS, Sₙ; colorrange=S_range, colormap=:haline) #note that this is still using old grid from T, S, initial, may need to recompute x and z using specific nodes 
Colorbar(fig[3, 2], hm_S, label="ppt")

# xρ, yρ, zρ = nodes(ρ_t[1])
# ax_ρ = Axis(fig[4, 1]; title="potential density", axis_kwargs...)
# hm_ρ = heatmap!(ax_ρ, xρ, zρ, ρₙ; colorrange=ρ_range, colormap=:viridis) #note that this is still using old grid from T, S, initial, may need to recompute x and z using specific nodes 
# Colorbar(fig[4, 2], hm_ρ, label="kg/m^3")

fig
@info "Making properties animation from data"
frames = 1:lastFrame
record(fig, filename * "properties.mp4", frames, framerate=8) do i
    @info string("Plotting frame ", i, " of ", frames[end])
    n[] = i
end





































# xT, yT, zT = nodes(T_t[1])
# ax_T = Axis(fig[3, 1]; title = "temperature", axis_kwargs...)
# hm_T = heatmap!(ax_T, xT, zT, Tₙ; colorrange = T_range, colormap = :thermal)
# Colorbar(fig[3,2], hm_T, label = "C")




#TODO: visualize streamlines & how they change

#visualize distributions
density_field = Field(density_operation)
compute!(density_field)

fig = Figure()
title = "initial conditions"
Label(fig[0, :], title)
x = xnodes(model.tracers.T)
z = znodes(model.tracers.T)

axis_kwargs = (xlabel="x (m)", ylabel="z (m)", width=400)

ax1 = Axis(fig[1, 1]; title="Temperature (C)", axis_kwargs...)
T_init_matrix = interior(model.tracers.T, :, 1, :)
T_lims = (minimum(T_init_matrix), maximum(T_init_matrix));
hm1 = heatmap!(ax1, x, z, T_init_matrix, colorrange=T_lims, colormap=:thermal, interpolate=true)
Colorbar(fig[1, 2], hm1)

ax2 = Axis(fig[2, 1]; title="Salinity (ppt)", axis_kwargs...)
S_init_matrix = interior(model.tracers.S, :, 1, :)
S_lims = (minimum(S_init_matrix), maximum(S_init_matrix));
hm2 = heatmap!(ax2, x, z, S_init_matrix, colorrange=S_lims, colormap=:haline, interpolate=true)
Colorbar(fig[2, 2], hm2)

ax3 = Axis(fig[3, 1]; title="Density (kg/m^3)", axis_kwargs...)
ρ_init_matrix = interior(density_field, :, 1, :)
ρ_lims = (minimum(ρ_init_matrix), maximum(ρ_init_matrix));
hm3 = heatmap!(ax3, x, z, ρ_init_matrix, colorrange=ρ_lims, colormap=Reverse(:viridis), interpolate=true)
Colorbar(fig[3, 2], hm3)
fig

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
x_spacings = xspacings(domain_grid, Center(), Face(), Center())
z_spacings = zspacings(domain_grid, Center(), Face(), Center())
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

x_spacings = xspacings(domain_grid, Center(), Face(), Center())
z_spacings = zspacings(domain_grid, Center(), Face(), Center())
