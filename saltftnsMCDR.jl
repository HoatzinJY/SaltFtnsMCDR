using Pkg
using Oceananigans
using OceanBioME
using TimeDates


#resolution parameters TODO determine best ones
#note: be careful of time resolution needed, spacing resoslution needed
#user set
grid_spacing = 0.1;# in meters
domain_x = 15; # width, in meters
domain_z = 15; #height, in meters

#calculated
x_res = domain_x/grid_spacing;
z_res = domain_z/grid_spacing;

#setting up model components

grid = RectilinearGrid(CPU(), Float64; size = (x_res, z_res), x = (0, domain_x), z = (-domain_z, 0), topology = (Bounded, Flat, Bounded))
clock = Clock{eltype(grid)}(time = 0);
advection = CenteredSecondOrder(); #default, not sure which one to choose
buoyancy = SeawaterBuoyancy(equation_of_state = LinearEquationofstate(thermal_expansion = 2e-4, haline_contraction = 78e-5)); #TODO: potentially add more accuracy here, currently set to global average
tracers = (:T, :S); #temperature, salinity
boundary_conditions::myBoundaries = myBoundaries();
timestepper = :QuasiAdamsBashforth2; #default, 3rd order option available
#not yet incorporated 
forcing::wallForcers = wallForcers(); #for imaginary pipe walls
immersed_boundary = nothing; #another option for imaginary pipe walls
closure = nothing; #for turbulent dissapation at edges of domain, can set to be direction and tracer specific, if diffusivity is nonconstant and calculated per timestep, need diffusivity_field as well
boundary_conditions::myBoundaries = myBoundaries(); # currently set to default (0 flux), but potentially could use immersed boundaries for
#biogeochemistry =  LOBSTER(; grid); #not yet used at all
background_fields::water_column_profiles = water_column_profiles(); #TODO, add in background t & s profiles to serve as a "contant" beyond pipe
#could add in hydrostatic_pressure_anomoly field to see if hydrostatic assumption is valid 

#sets up model
#following are considered negligable/not accounted for: coriolis, stokes drift
#following are determinined automatically: pressure solver 
#have yet to think about following: auxiliary fields
model = NonhydrostaticModel(; grid, clock, advection, buoyancy, boundary_conditions, tracers, timestepper)

#setting initial conditions 
#TODO: perturb one part, maybe give it an intial velocity from pumping?? (this would align with matlab)
