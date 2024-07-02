using Pkg
using Oceananigans
using OceanBioME
using TimesDates
using CairoMakie
using Printf
using Oceananigans, SeawaterPolynomials.TEOS10
using Oceananigans.Models: seawater_density

const day = 86400;
const hour = 3600; 
const minute = 60;
geopotential_height = 0; # sea surface height for potential density calculations

#pipe parameters
pipe_radius = 0.5;
pipe_length = 10;
pipe_top_depth = 1;

#pumping parameters
height_displaced = 3;
initial_pipe_velocity = 0.001; #not yet used 

#resolution parameters TODO determine best ones
#note: be careful of time resolution needed, spacing resoslution needed
#user set
grid_spacing = 0.1;# in meters
domain_x = 15; # width, in meters
domain_z = 15; #height, in meters

#calculated
x_res = floor(Int, domain_x/grid_spacing);
z_res = floor(Int, domain_z/grid_spacing);
locs = (Center(), Face(), Center());

#setting up model components
domain_grid = RectilinearGrid(CPU(), Float64; size = (x_res, z_res), x = (0, domain_x), z = (-domain_z, 0), topology = (Bounded, Flat, Bounded))
clock = Clock{eltype(domain_grid)}(time = 0);
advection = CenteredSecondOrder(); #default, not sure which one to choose
eos = TEOS10EquationOfState()
buoyancy = SeawaterBuoyancy(equation_of_state = eos)
#buoyancy = SeawaterBuoyancy(equation_of_state = LinearEquationOfState(thermal_expansion = 2e-4, haline_contraction = 78e-5)); #TODO: potentially add more accuracy here, currently set to global average
tracers = (:T, :S); #temperature, salinity
#ScalarDiffusivity(ν=1e-6, κ=(S=1e-7, T=1.45e-10)) #this throws a negative under sq rt error in TEOS-10 TODO; set this/make it work, get more accurate constants from https://web.mit.edu/seawater/2017_MIT_Seawater_Property_Tables_r2b_2023c.pdf
timestepper = :QuasiAdamsBashforth2; #default, 3rd order option available
#not yet incorporated 
# forcing = wallForcers(); #for imaginary pipe walls
# closure = nothing; #for turbulent dissapation at edges of domain, can set to be direction and tracer specific, if diffusivity is nonconstant and calculated per timestep, need diffusivity_field as well
# #biogeochemistry =  LOBSTER(; domain_grid); #not yet used at all
# background_fields::water_column_profiles = water_column_profiles(); #TODO, add in background t & s profiles to serve as a "contant" beyond pipe
# #could add in hydrostatic_pressure_anomoly field to see if hydrostatic assumption is valid 

#sets up model
#following are considered negligable/not accounted for: coriolis, stokes drift
#following are determinined automatically: pressure solver 
#have yet to think about following: auxiliary fields

model = NonhydrostaticModel(; grid = domain_grid, clock, advection, buoyancy, tracers, timestepper)

#helper functions
density_operation = seawater_density(model; geopotential_height)

#setting up just a basic gradient for water column 
T_top = 21.67;
T_bot = 11.86;
S_bot = 35.22;
S_top = 34.18;
delta_z = 200;
#sets initial t and s 
x_center = domain_x/2;
function T_init(x, z)
    #inside pipe
    if (x > (x_center - pipe_radius) && x < (x_center + pipe_radius) && z < -pipe_top_depth && z > -(pipe_top_depth + pipe_length))
        return T_top - ((T_bot - T_top)/delta_z)*(z + height_displaced); 
    #outside pipe
    else 
        return T_top - ((T_bot - T_top)/delta_z)z;
    end
end
function S_init(x, z)
    if (x > (x_center - pipe_radius) && x < (x_center + pipe_radius) && z < -pipe_top_depth && z > -(pipe_top_depth + pipe_length))
        return S_top - ((S_bot - S_top)/delta_z)*(z + height_displaced); 
    else 
        return S_top - ((S_bot - S_top)/delta_z)z;
    end
end
function w_init(x, z)
    if (x > (x_center - pipe_radius) && x < (x_center + pipe_radius) && z < -pipe_top_depth && z > -(pipe_top_depth + pipe_length))
        #return initial_pipe_velocity;
        return 0; 
    else 
        return 0;
    end
end 
set!(model, T = T_init, S = S_init, w = w_init)

#running model
#TODO: figure out why it is throwing a complext domain error, line 83, when calculating salinity in TEOS-10
simulation = Simulation(model, Δt = 1, stop_time = 15day, wall_time_limit = 180) # make initial delta t bigger
timeWizard = TimeStepWizard(cfl = 0.33, max_Δt = 10minute) #TODO: set max delta t?
simulation.callbacks[:timeWizard] = Callback(timeWizard, IterationInterval(4)) 
#log progress --> TODO: add this

#fields = Dict("u" => model.velocities.u, "w" => model.velocities.w, "T" => model.tracers.T, "S" => model.tracers.S)
w = model.velocities.w;
u = model.velocities.u;
T = model.tracers.T;
S = model.tracers.S;
ζ = Field( - ∂x(w) + ∂z(u)) #vorticity in y 
ρ = Field(density_operation)
filename = "tryseven"
simulation.output_writers[:outputs] = JLD2OutputWriter(model, (; u, w, T, S, ζ, ρ); filename, schedule = TimeInterval(10), overwrite_existing = true)
#time average perhaps?

run!(simulation; pickup = false)

#visualize simulation
#may want to average

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
title = @lift @sprintf("t = %1.2f days", round(times[$n] / day, digits=2))
uₙ = @lift interior(u_t[$n], :, 1, :)
wₙ = @lift interior(w_t[$n], :, 1, :)
ζₙ = @lift interior(ζ_t[$n], :, 1, :)
Tₙ = @lift interior(T_t[$n], :, 1, :)
Sₙ = @lift interior(S_t[$n], :, 1, :)
ρₙ = @lift interior(ρ_t[$n], :, 1, :)

num_Data_Points = size(u_t.data, 4)
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
ρ_range = getMaxAndMin(num_Data_Points, ρ_t)
ζ_range = getMaxAndMin(num_Data_Points, ζ_t)


fig = Figure(size = (600, 900))
axis_kwargs = (xlabel = "x (m)", ylabel = "z (m)", width = 400)
fig[1, :] = Label(fig, title)

xw, yw, zw = nodes(w_t[1])
ax_w = Axis(fig[2, 1]; title = "w velocity", axis_kwargs...)
hm_w = heatmap!(ax_w, xw, zw, wₙ; colorrange = w_range, colormap = :balance) #note that this is still using old grid from T, S, initial, may need to recompute x and z using specific nodes 
Colorbar(fig[2,2], hm_w, label = "m/s")

xu, yu, zu = nodes(u_t[1])
ax_u = Axis(fig[3, 1]; title = "u velocity", axis_kwargs...)
hm_u = heatmap!(ax_u, xu, zu, uₙ; colorrange = u_range, colormap = :balance) #note that this is still using old grid from T, S, initial, may need to recompute x and z using specific nodes 
Colorbar(fig[3,2], hm_u, label = "m/s")

xζ, yζ, zζ = nodes(ζ_t[1])
ax_ζ = Axis(fig[4, 1]; title = "vorticity", axis_kwargs...)
hm_ζ = heatmap!(ax_ζ, xζ, zζ, ζₙ; colorrange = ζ_range, colormap = :balance) #note that this is still using old grid from T, S, initial, may need to recompute x and z using specific nodes 
Colorbar(fig[4,2], hm_ζ, label = "rot/sec")

fig
# xT, yT, zT = nodes(T_t[1])
# ax_T = Axis(fig[3, 1]; title = "temperature", axis_kwargs...)
# hm_T = heatmap!(ax_T, xT, zT, Tₙ; colorrange = T_range, colormap = :thermal)
# Colorbar(fig[3,2], hm_T, label = "C")



@info "Making animation from data"
frames = 1:num_Data_Points #fix to length[]
record(fig, filename * ".mp4", frames, framerate = 16) do i 
    @info string("Plotting frame ", i, " of ", frames[end])
    n[] = i
end

#TODO: visualize streamlines & how they change

#visualize distributions
density_field = Field(density_operation)
compute!(density_field)

fig = Figure()
title = "initial conditions"
Label(fig[0, :], title)
x = xnodes(model.tracers.T)
z = znodes(model.tracers.T)

axis_kwargs = (xlabel = "x (m)", ylabel = "z (m)", width = 400)

ax1 = Axis(fig[1, 1]; title = "Temperature (C)", axis_kwargs...)
T_init_matrix = interior(model.tracers.T, :, 1, :)
T_lims = (minimum(T_init_matrix), maximum(T_init_matrix));
hm1 = heatmap!(ax1, x , z, T_init_matrix, colorrange = T_lims, colormap = :thermal, interpolate = true)
Colorbar(fig[1, 2], hm1)

ax2 = Axis(fig[2, 1]; title = "Salinity (ppt)", axis_kwargs...)
S_init_matrix = interior(model.tracers.S, :, 1, :)
S_lims = (minimum(S_init_matrix), maximum(S_init_matrix));
hm2 = heatmap!(ax2, x , z, S_init_matrix, colorrange = S_lims, colormap = :haline, interpolate = true)
Colorbar(fig[2, 2], hm2)

ax3 = Axis(fig[3, 1]; title = "Density (kg/m^3)", axis_kwargs...)
ρ_init_matrix = interior(density_field, :, 1, :)
ρ_lims = (minimum(ρ_init_matrix), maximum(ρ_init_matrix));
hm3 = heatmap!(ax3, x , z, ρ_init_matrix, colorrange = ρ_lims, colormap = Reverse(:viridis), interpolate = true)
Colorbar(fig[3, 2], hm3)
fig

#check that model has reset - plot velocities
fig = Figure()
title = "initial velocities"
Label(fig[0, :], title)
x = xnodes(model.tracers.T)
z = znodes(model.tracers.T)

ax1 = Axis(fig[2, 1]; title = "u velocity (m/s)", axis_kwargs...)
u_init_matrix = interior(model.velocities.u, :, 1, :)
u_lims = (minimum(u_init_matrix), maximum(u_init_matrix));
hm1 = heatmap!(ax1, x , z, u_init_matrix, colorrange = u_lims, colormap = Reverse(:balance), interpolate = true)
Colorbar(fig[2, 2], hm1)

axis_kwargs = (xlabel = "x (m)", ylabel = "z (m)", width = 400)
ax2 = Axis(fig[1, 1]; title = "w velocity (m/s)", axis_kwargs...)
w_init_matrix = interior(model.velocities.w, :, 1, :)
w_lims = (minimum(w_init_matrix), maximum(w_init_matrix));
hm2 = heatmap!(ax2, x , z, w_init_matrix, colorrange = w_lims, colormap = Reverse(:balance), interpolate = true)
Colorbar(fig[1, 2], hm2)

display(fig)




#testing
#setting up pipe
# domain_arr = nodes(domain_grid, locs; reshape = true)
# pipe_interior_grid_x = AbstractArray{Float64, 3}

# x_spacings = xspacings(domain_grid, Center(), Face(), Center())
# z_spacings = zspacings(domain_grid, Center(), Face(), Center())
