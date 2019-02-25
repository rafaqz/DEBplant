using Revise, Unitful, Microclimate, JLD2, DataStructures
using Unitful: W, m, °C, hr, mol, K, s, J, Mg, kg, kPa, Pa

basepath = ENV["MICROCLIM"]
shade = 0
years = 2001:2011

skipped = (:snowdepth, :soilwatercontent) 
envgrid = load_grid(basepath, years, shade, skipped);
# tas = MicroclimPoint(envgrid, CartesianIndex(56, 53))
# desert = MicroclimPoint(envgrid, CartesianIndex(30,25))
# qld = MicroclimPoint(envgrid, CartesianIndex(66, 29))

t1 = MicroclimPoint(envgrid, CartesianIndex(65, 35))
t2 = MicroclimPoint(envgrid, CartesianIndex(60, 35))
t3 = MicroclimPoint(envgrid, CartesianIndex(55, 35))
t4 = MicroclimPoint(envgrid, CartesianIndex(50, 35))

environments = Dict(:t1 => t1, :t2 => t2, :t3 => t3, :t4 => t4)

@save "microclimate/locations.jld" environments
# @load "microclimate/locations.jld" tas desert qld
