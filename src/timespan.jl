using Unitful, Microclimate, JLD2
using Unitful: W, m, °C, hr, mol, K, s, J, Mg, kg, kPa, Pa

basepath = ENV["MICROCLIM"]

environment = load_point(basepath, 1995:1996, 0, CartesianIndex(5, 20))

@save "timespan.jld" environment
