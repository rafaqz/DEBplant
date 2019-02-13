using Unitful, Microclimate, JLD2
using Unitful: W, m, °C, hr, mol, K, s, J, Mg, kg, kPa, Pa

basepath = ENV["MICROCLIM"]
shade = 0
years = 2001:2002

envgrid = load_grid(basepath, years, shade)
tas = MicroclimPoint(envgrid, CartesianIndex(56, 53))
desert = MicroclimPoint(envgrid, CartesianIndex(30,25))
qld = MicroclimPoint(envgrid, CartesianIndex(66, 29))

@save "locations.jld" tas desert qld
@load "locations.jld" tas desert qld
