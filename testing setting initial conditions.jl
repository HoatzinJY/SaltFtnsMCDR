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
using NCDatasets
using Dierckx

scale = 1
x_grid_spacing = 0.05scale;# in meters
z_grid_spacing = 0.05scale;# in meters
domain_x = 15scale; # width, in meters
domain_z = 15scale; #height, in meters
roundUp(num::Float64, base) = ceil(Int, num / base) * base\
findNearest(A::AbstractArray, x) = findmin(abs(A-x)) 

#pipe parameters
pipe_radius = 0.5scale;
pipe_length = 10scale;
pipe_top_depth = 1scale;
pipe_wall_thickness_intended = 0.01scale;
pipe_wall_thickness = roundUp(pipe_wall_thickness_intended, x_grid_spacing)
@info @sprintf("Pipe walls are %1.2f meters thick", pipe_wall_thickness)

#setting initial conditions
#TODO: be able to set filtname by month
#specify month, latitude and longitude of interest 
month = "06"
lat_p = 22
lon_p = -150 

#open files
T_filename = joinpath("WOA18", "woa18_decav_t"*month*"_01.nc")
S_filename = joinpath("WOA18", "woa18_decav_s"*month*"_01.nc")
T_file = NetCDF.open(T_filename)

#find location that corresponds to lat/long
lat = NetCDF.readvar(T_file["lat"])
lon = NetCDF.readvar(T_file["lon"])
depth = NetCDF.readvar(T_file["depth"])
NetCDF.close(T_file)
lat_Index = findNearest(lat, lat_p)[1]
lon_Index = findNearest(lon, lon_p)[1]



function T_init(x, z)

end
T_init(x, y, z) = T_init(x, z)
function S_init(x, z)

end
S_init(x, y, z) = S_init(x, z)



domain_grid = RectilinearGrid(CPU(), Float64; size=(x_res, z_res), x=(0, domain_x), z=(-domain_z, 0), topology=(Bounded, Flat, Bounded))
clock = Clock{eltype(domain_grid)}(time=0);
advection = CenteredSecondOrder()
buoyancy = nothing;
#buoyancy = SeawaterBuoyancy(equation_of_state = LinearEquationOfState(thermal_expansion = 2e-4, haline_contraction = 78e-5)); #TODO: potentially add more accuracy here, currently set to global average
tracers = (:T, :S); #temperature, salinity
timestepper = :QuasiAdamsBashforth2; #default, 3rd order option available 
closure = ScalarDiffusivity(ν=viscosity.molecular, κ=(T=T_diffusivity.molecular, S=S_diffusivity.molecular))
noforcing(x, z, t) = 0 
forcing = (u = noforcing,  w = noforcing)
model = NonhydrostaticModel(; grid=domain_grid, clock, advection, buoyancy, tracers, timestepper, closure, forcing)
@info "model made"
#helper functions
density_operation = seawater_density(model; geopotential_height)

set!(model, T=T_init, S=S_init, w=w_init)
@info "initial ocnditions set"

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