using Unitful
using Unitful: W, m, °C, hr, mol, K, s, J, Mg, kg, kPa, Pa
using DynamicEnergyBudgets: Env
using NCDatasets
using JLD2

basepath = "/home/raf/Uni/Masters/NetCDF/"

load_netcdf(fpath) = Dataset(joinpath(basepath, fpath))

loadall(filename, varname, index) = begin
    println("adding ", varname)
    all = Float64[]
    for year in 1990:2017
        annual = load_netcdf(filename * string(year) * ".nc")[varname][index, :]
        all = vcat(all, annual)
    end
    all
end

build_env(rad, tair, wind, relh, tsoil, pot, i::CartesianIndex) = begin
    Env(rad .* 0.1 * W*m^-2, 
        map(x -> x ./ 10 .* °C .|> K, tair), 
        map(x -> x ./ 100, relh), 
        map(x -> x ./ 10 .* m*s^-1, wind), 
        map(x -> x ./ 10 * °C .|> K, tsoil),  
        map(x -> x ./ 10 * J * kg^-1 * 1Mg * m^-3 .|> kPa, pot)
       )  
end


solr = loadall("SOLR/SOLR_", "SOLR", i)
rh1cm = loadall("RH1cm_0pctShade/RH1cm_0pctShade_", "RH1cm", i)
rh120cm = loadall("RH120cm/RH120cm_", "RH120cm", i)
ta1cm = loadall("Tair1cm_0pctShade/TA1cm_0pctShade_", "TA1cm", i)
ta120cm = loadall("Tair120cm/TA120cm_", "TA120cm", i)
v1cm = loadall("V1cm/V1cm_", "V1cm", i)
v120cm = loadall("V120cm/V120cm_", "V120cm", i)
pot0cm = loadall("pot0cm_0pctShade/pot0cm_0pctShade_", "pot0cm", i)
pot2_5cm = loadall("pot2.5cm_0pctShade/pot2.5cm_0pctShade_", "pot2.5cm", i)
pot5cm = loadall("pot5cm_0pctShade/pot5cm_0pctShade_", "pot5cm", i)
pot10cm = loadall("pot10cm_0pctShade/pot10cm_0pctShade_", "pot10cm", i)
pot20cm = loadall("pot20cm_0pctShade/pot20cm_0pctShade_", "pot20cm", i)
pot30cm = loadall("pot30cm_0pctShade/pot30cm_0pctShade_", "pot30cm", i)
pot50cm = loadall("pot50cm_0pctShade/pot50cm_0pctShade_", "pot50cm", i)
pot100cm = loadall("pot100cm_0pctShade/pot100cm_0pctShade_", "pot100cm", i)
soil0cm = loadall("soil5cm_0pctShade/soil5cm_0pctShade_", "soil5cm", i)
soil2_5cm = loadall("soil2.5cm_0pctShade/soil2.5cm_0pctShade_", "soil2.5cm", i)
soil5cm = loadall("soil0cm_0pctShade/soil0cm_0pctShade_", "soil0cm", i)
soil10cm = loadall("soil10cm_0pctShade/soil10cm_0pctShade_", "soil10cm", i)
soil20cm = loadall("soil20cm_0pctShade/soil20cm_0pctShade_", "soil20cm", i)
soil30cm = loadall("soil30cm_0pctShade/soil30cm_0pctShade_", "soil30cm", i)
soil50cm = loadall("soil50cm_0pctShade/soil50cm_0pctShade_", "soil50cm", i)
soil100cm = loadall("soil100cm_0pctShade/soil100cm_0pctShade_", "soil100cm", i)

rad = solr
tair = ta1cm, ta120cm
wind = v1cm, v120cm
relh = rh1cm, rh120cm
tsoil = soil0cm, soil2_5cm, soil5cm, soil10cm, soil20cm, soil30cm, soil50cm, soil100cm
pot = pot0cm, pot2_5cm, pot5cm, pot10cm, pot20cm, pot30cm, pot50cm, pot100cm

environment = build_env(rad, tair, wind, relh, tsoil, pot, CartesianIndex(5, 20))

i = CartesianIndex(20, 30)

@save "long_env.jld" environment
