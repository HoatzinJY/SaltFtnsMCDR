using Pkg
using Oceananigans
using OceanBioME
using TimesDates
using CairoMakie

#pipe parameters
pipe_radius = 0.5;
pipe_length = 10;
pipe_top_depth = 1;

#resolution parameters TODO determine best ones
#note: be careful of time resolution needed, spacing resoslution needed
#user set
grid_spacing = 0.1;# in meters
domain_x = 15; # width, in meters
domain_z = 15; #height, in meters

#calculated
x_res = floor(Int, domain_x/grid_spacing);
z_res = floor(Int, domain_z/grid_spacing);

#setting up model components
domain_grid = RectilinearGrid(CPU(), Float64; size = (x_res, z_res), x = (0, domain_x), z = (-domain_z, 0), topology = (Bounded, Flat, Bounded))
pipe_grid = RectilinearGrid();
clock = Clock{eltype(domain_grid)}(time = 0);
advection = CenteredSecondOrder(); #default, not sure which one to choose
buoyancy = SeawaterBuoyancy(equation_of_state = LinearEquationOfState(thermal_expansion = 2e-4, haline_contraction = 78e-5)); #TODO: potentially add more accuracy here, currently set to global average
tracers = (:T, :S); #temperature, salinity
timestepper = :QuasiAdamsBashforth2; #default, 3rd order option available
#not yet incorporated 
# forcing = wallForcers(); #for imaginary pipe walls
# immersed_boundary = nothing; #another option for imaginary pipe walls
# closure = nothing; #for turbulent dissapation at edges of domain, can set to be direction and tracer specific, if diffusivity is nonconstant and calculated per timestep, need diffusivity_field as well
# boundary_conditions= myBoundaries(); # currently set to default (0 flux), but potentially could use immersed boundaries for
# #biogeochemistry =  LOBSTER(; domain_grid); #not yet used at all
# background_fields::water_column_profiles = water_column_profiles(); #TODO, add in background t & s profiles to serve as a "contant" beyond pipe
# #could add in hydrostatic_pressure_anomoly field to see if hydrostatic assumption is valid 

#sets up model
#following are considered negligable/not accounted for: coriolis, stokes drift
#following are determinined automatically: pressure solver 
#have yet to think about following: auxiliary fields
model = NonhydrostaticModel(; domain_grid, clock, advection, buoyancy, tracers, timestepper)

#setting up just a basic gradient for water column 
T_top = 21.67;
T_bot = 11.86;
S_bot = 35.22;
S_top = 34.18;
delta_z = 200;
#sets up initial gradient, starting from surface
#TODO: pipe, can use for loop to set up mask, or somehow find a mask
T_initial(x, z) = T_top - ((T_bot - T_top)/delta_z)z;
S_initial(x, z) = S_top - ((S_bot - S_top)/delta_z)z;
set!(model, T = T_initial; S = S_initial)


#visualize inital temperature distribution
fig = Figure()
title = "initial conditions"
Label(fig[0, :], title)

axis_kwargs = (xlabel = "x (m)", ylabel = "z (m)", width = 400)
ax1 = Axis(fig[1, 1]; title = "Temperature (C)", axis_kwargs...)
ax2 = Axis(fig[2, 1]; title = "Salinity (ppt)", axis_kwargs...)

T_lims = (T_initial(0, -domain_z), T_top);
S_lims = (S_top, S_initial(0, -domain_z));
x = xnodes(model.tracers.T)
z = znodes(model.tracers.T)
T_init_matrix = interior(model.tracers.T, :, 1, :)
S_init_matrix = interior(model.tracers.S, :, 1, :)
hm1 = heatmap!(ax1, x , z, T_init_matrix, colorrange = T_lims, colormap = :thermal, interpolate = true)
hm2 = heatmap!(ax2, x , z, S_init_matrix, colorrange = S_lims, colormap = :haline, interpolate = true)
Colorbar(fig[1, 2], hm1)
Colorbar(fig[2, 2], hm2)

fig 

#visualize initial salinity distribution


#setting initial conditions 
#TODO: perturb one part, maybe give it an intial velocity from pumping?? (this would align with matlab)