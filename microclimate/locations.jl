# Build location microclimates from Microclim datasets

using Revise, Unitful, Microclimate, JLD2, DataStructures
using Unitful: W, m, °C, hr, mol, K, s, J, Mg, kg, kPa, Pa

basepath = ENV["MICROCLIM"]
shade = 0
years = 2005:2011

skipped = (:soilwatercontent,) 
envgrid = load_grid(basepath, years, shade, skipped);
# tas = MicroclimPoint(envgrid, CartesianIndex(56, 53))
# desert = MicroclimPoint(envgrid, CartesianIndex(30,25))
# qld = MicroclimPoint(envgrid, CartesianIndex(66, 29))

t1 = MicroclimPoint(envgrid, CartesianIndex(65, 35))
t2 = MicroclimPoint(envgrid, CartesianIndex(60, 35))
t3 = MicroclimPoint(envgrid, CartesianIndex(55, 35))

@save "microclimate/locations.jld" t1 t2 t3
# @load "microclimate/locations.jld" tas desert qld
