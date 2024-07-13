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

"""currently only uses molecular diffusion as there is no motion""";

#= 
NOTES
- seems to be significanty slower after adding in forcing 
PROBLEM: not setting diffusive cfl
- ok at 0.2, bad at 0.33 *diffusive time scale 
- TODO: figure ouhow to use difusive cfl
=#

#TODO: test that velocity forcing has no impact on model 
#TODO: test with buoyancy
#TODO: forcing: no hoizontal or vertical eddy diffusion through pipe walls, only molecular 
const day = 86400;
const hour = 3600;
const minute = 60;
const g = 9.806; #gravitational acceleration, in m^2/s
myCFL = 0.33
scale = 1

geopotential_height = 0; # sea surface height for potential density calculations
roundUp(num::Float64, base) = ceil(Int, num / base) * base

T_top = 21.67;
T_bot = 11.86;
S_bot = 35.22;
S_top = 34.18;
delta_z = 200;

#grid spacing
x_grid_spacing = 0.05scale;# in meters
z_grid_spacing = 0.05scale;# in meters
domain_x = 15scale; # width, in meters
domain_z = 15scale; #height, in meters

#pipe parameters
pipe_radius = 0.5scale;
pipe_length = 10scale;
pipe_top_depth = 1scale;
pipe_wall_thickness_intended = 0.01scale;
pipe_wall_thickness = roundUp(pipe_wall_thickness_intended, x_grid_spacing)
@info @sprintf("Pipe walls are %1.2f meters thick", pipe_wall_thickness)

#pumping parameters
height_displaced = 3scale;
initial_pipe_velocity = 0.001scale; #not yet used 

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
function  waterBorderMask(x, z)
    if(min(x, domain_x - x) < x_grid_spacing || (domain_z + z) < z_grid_spacing)
        return 1
    else
        return 0
    end
end
waterBorderMask(x, y, z) = waterBorderMask(x, z)


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
T_init(x, y, z) = T_init(x, z)
function S_init(x, z)
    if (isInsidePipe(x, z))
        return S_top - ((S_bot - S_top) / delta_z) * (z + height_displaced)
    else
        return S_top - ((S_bot - S_top) / delta_z)z
    end
end
S_init(x, y, z) = S_init(x, z)

domain_grid = RectilinearGrid(CPU(), Float64; size=(x_res, z_res), x=(0, domain_x), z=(-domain_z, 0), topology=(Bounded, Flat, Bounded))
advection = CenteredSecondOrder();
tracers = (:T, :S); 
timestepper = :RungeKutta3; 


function tempDiffusivities(x, y, z, parameters::NamedTuple, wall_indicator::String)
    if (z < 9 || wall_indicator == "WALL")
        return 1 #edit to account for wall thickness
    else 
        return 0
    end
end
#to test, with pipe wall, this sets top half pipe wall with big diffusivity, rest 2 degrees of magnitude less
# function tempDiffusivities(x, y, z, parameters::NamedTuple, wall_indicator::String)
#     if ((isPipeWall(x, z) && z < 9) || wall_indicator == "WALL")
#         return 1
#     else 
#         return 0.01
#     end
# end
tempDiffusivities(x, z) = tempDiffusivities(x, 0, z, diffusivity_data, "")
tempDiffusivities(x, y, z) = tempDiffusivities(x, y, z, diffusivity_data, "")
tempDiffusivities(x, y, z, parameters::NamedTuple) = tempDiffusivities(x, y, z, parameters, "")
#option for walls - currently only thermal
closure = ScalarDiffusivity(ν=viscosity.molecular, κ=(T=tempDiffusivities, S=1.3E-9))

#closure = ScalarDiffusivity(ν=0, κ=(T=1.46e-7, S=1.3E-9))
#closure = ScalarDiffusivity(ν=0, κ=(S=0.5, T=1.0))


buoyancy = nothing
# u_damping_rate = 1/0.1 #relaxes fields on 0.1 second time scale
# w_damping_rate = 1/1 #relaxes fields on 1 second time scale
# u_pipe_wall = Relaxation(rate = u_damping_rate, mask = pipeWallMask)
# w_pipe_wall = Relaxation(rate = w_damping_rate, mask = pipeWallMask)
# forcing = (u = (u_pipe_wall),  w = (w_pipe_wall))

noforcing(x, y, z) = 0 
forcing = (u = noforcing,  w = noforcing)
# wall_damping_rate = 1/0.000001
# T_border = Relaxation(rate = wall_damping_rate, mask = waterBorderMask, target = T_init)
# S_border = Relaxation(rate = wall_damping_rate, mask = waterBorderMask, target = S_init)
# forcing = (T = T_border, S=S_border)

model = NonhydrostaticModel(; grid = domain_grid, advection, buoyancy, tracers, timestepper, closure, forcing)
@info "model made"
set!(model, T=T_init, S=S_init)
@info "initial ocnditions set"

min_grid_spacing = min(minimum_xspacing(model.grid), minimum_zspacing(model.grid))
diffusion_time_scale = (min_grid_spacing^2)/model.closure.κ.T
initial_time_step = 0.1*diffusion_time_scale
#max_time_step = myCFL*diffusion_time_scale
simulation_duration = 15day
run_duration = 2minute;

simulation = Simulation(model, Δt=initial_time_step, stop_time=simulation_duration, wall_time_limit=run_duration) # make initial delta t bigger
timeWizard = TimeStepWizard(cfl=0.2, diffusive_cfl = 0.2) #TODO: set max delta t?
simulation.callbacks[:timeWizard] = Callback(timeWizard, IterationInterval(4))
progress_message(sim) = @printf("Iteration: %04d, time: %s, Δt: %s, wall time: %s\n", iteration(sim), prettytime(sim), prettytime(sim.Δt), prettytime(sim.run_wall_time))
add_callback!(simulation, progress_message, IterationInterval(10))


T = model.tracers.T;
S = model.tracers.S;
filename = "diffusion trial with cfldiffusion on to 0.2"
simulation.output_writers[:outputs] = JLD2OutputWriter(model, (; T, S); filename, schedule= IterationInterval(10), overwrite_existing=true) #can also set to TimeInterval
#time average perhaps?

run!(simulation; pickup=false)

output_filename = filename * ".jld2"

T_t = FieldTimeSeries(output_filename, "T")
S_t = FieldTimeSeries(output_filename, "S")
times = T_t.times
n = Observable(1)
Observable(1)
Tₙ = @lift interior(T_t[$n], :, 1, :)
Sₙ = @lift interior(S_t[$n], :, 1, :)
num_Data_Points = length(times)
num_Data_Points = 4
function getMaxAndMin(numPoints, dataSeries)
    myMax = maximum(interior(dataSeries[1], :, 1, :))
    myMin = minimum(interior(dataSeries[1], :, 1, :))
    for i in 2:numPoints
        myMax = max(maximum(interior(dataSeries[i], :, 1, :)),myMax )
        myMin = min(minimum(interior(dataSeries[i], :, 1, :)), myMin)
    end
    return (myMin, myMax)
end
T_range = getMaxAndMin(num_Data_Points, T_t)
S_range = getMaxAndMin(num_Data_Points, S_t)


#animate properties
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

fig
@info "Making properties animation from data"
frames = 1:length(times)
record(fig, filename * "properties.mp4", frames, framerate=8) do i
    @info string("Plotting frame ", i, " of ", frames[end])
    n[] = i
end