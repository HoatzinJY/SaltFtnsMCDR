using Pkg
using Oceananigans
using OceanBioME
using TimesDates
using CairoMakie
using Printf
using SeawaterPolynomials
using GibbsSeaWater
using NetCDF

#plotting functions
function plotInitialConditions(m)
    set_theme!(Theme(fontsize = 24, linewidth=3))

    fig = Figure()
    title = "initial Conditions"
    Label(fig[0, :], title)

    x = xnodes(model.tracers.A)
    A = interior(model.tracers.A, :, 1, 1)

    axis = Axis(fig[2, 1]; xlabel = "x(m)", ylabel = "Tracer(C)", width = 400)
    lines1 = lines!(axis, x, A)
    fig
end 

#FIGURE OUT WHAT A DIFFUSIVITY FIEDLD IS? 

function maskCenter(x)
    if (1 < x < 2)
        return 1
    else
        return 0
    end
end

function A_init(x)
    return 5
end    

#THINGS TO DEFINE  


trial_name = "callbacktestone"

grid = RectilinearGrid(size=3, x=(0, 3), topology=(Bounded, Flat, Flat))
clock = Clock{eltype(grid)}(time=0)
timestepper = :QuasiAdamsBashforth2;
tracers =:A
closure = ScalarDiffusivity(ν=0, κ = 1.0)

A_bcs = FieldBoundaryConditions(east = ValueBoundaryCondition(10), west = ValueBoundaryCondition(1))
boundary_conditions = (A = A_bcs,)

function A_forcing_func(x, t, A, time_step)
    @info @sprintf("In function, x = %.2f | A = %.2f | ΔT = %s", x, A, prettytime(time_step))
    return maskCenter(x) * 1/(1.01*time_step) * (0 - A)
end 

diffusion_time_scale = (1^2)/1
initial_time_step = 0.1*diffusion_time_scale 
time_step = initial_time_step
A_forcing = Forcing(A_forcing_func, parameters=time_step, field_dependencies=:A)
forcing = (A = A_forcing,)

model = NonhydrostaticModel(; grid, clock, timestepper, tracers, closure, boundary_conditions, forcing)
set!(model, A = A_init)
@info "initial conditions set"
plotInitialConditions(model)


num_iterations = 3
run_duration = 30
simulation = Simulation(model, Δt = initial_time_step, stop_iteration = num_iterations, wall_time_limit = run_duration)
timeWizard = TimeStepWizard(cfl=0.2, diffusive_cfl = 0.2) 
# function timeWizardMessage()
#     @info @sprintf("*****WIZARD CALLED*****, Δt = %s", prettytime(sim.Δt))
#     timeWizard
#     @info @sprintf("*****WIZARD ENDED*****, Δt = %s", prettytime(sim.Δt))
# end
#TODO: work on gettng those to call time 
simulation.callbacks[:timeWizard] = Callback(timeWizard, IterationInterval(1))
progress_message_timestep(sim) = @printf("TIMESTEP Iteration: %04d, time: %.3e seconds, Δt: %s\n", iteration(sim), time(sim), prettytime(sim.Δt))
progress_message_tendency(sim) = @printf("TENDENCY Iteration: %04d, time: %.3e seconds, Δt: %s\n", iteration(sim), time(sim), prettytime(sim.Δt))
progress_message_updatestate(sim) = @printf("UPDATESTATE Iteration: %04d, time: %s, Δt: %s\n", iteration(sim), prettytime(sim), prettytime(sim.Δt))
add_callback!(simulation, progress_message_tendency, IterationInterval(1), callsite = TendencyCallsite())
add_callback!(simulation, progress_message_updatestate, IterationInterval(1), callsite = UpdateStateCallsite())
add_callback!(simulation, progress_message_timestep, IterationInterval(1), callsite = TimeStepCallsite())


# iteration(simulation)
# time(simulation)
# simulation.Δt
#changing model forcing


#saving data
A = model.tracers.A
filename = joinpath("callback_testing",trial_name)
simulation.output_writers[:outputs] = JLD2OutputWriter(model, (; A); filename, schedule=IterationInterval(1), overwrite_existing=true) #can also set to TimeInterval
run!(simulation; pickup = false)

#visualize
output_filename = filename * ".jld2"
A_t = FieldTimeSeries(output_filename, "A")
times = A_t.times
n = Observable(1)
Observable(1)
Aₙ = @lift interior(A_t[$n], :, 1, 1)
numPoints = length(times)

#plot 
fig = Figure(size = (400, 300))
title = @lift @sprintf("t = %s", prettytime(times[$n]))
axis_kwargs = (xlabel="x (m)", ylabel="Tracer", width=400)
fig[1, :] = Label(fig, title)
xA, yA, zA = nodes(A_t[1])
ax_A = Axis(fig[2, 1]; title="tracer vs distance", axis_kwargs...)
lines_A = lines!(ax_A, xA, Aₙ) #note that this is still using old grid from T, S, initial, may need to recompute x and z using specific nodes 
fig
@info "Making properties animation from data"
frames = 1:numPoints
record(fig, filename * "video.mp4", frames, framerate=1) do i
    @info string("Plotting frame ", i, " of ", frames[end])
    n[] = i
end