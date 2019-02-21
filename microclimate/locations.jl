using Revise, Unitful, Microclimate, JLD2
using Unitful: W, m, °C, hr, mol, K, s, J, Mg, kg, kPa, Pa

basepath = ENV["MICROCLIM"]
shade = 0
years = 2001:2011

skipped = (:snowdepth, :windspeed, :relhumidity, :soilwatercontent) 
envgrid = load_grid(basepath, years, shade, skipped);
tas = MicroclimPoint(envgrid, CartesianIndex(56, 53))
desert = MicroclimPoint(envgrid, CartesianIndex(30,25))
qld = MicroclimPoint(envgrid, CartesianIndex(66, 29))
soilwaterpotential(tas)

@save "microclimate/locations.jld" tas desert qld
# @load "microclimate/locations.jld" tas desert qld
