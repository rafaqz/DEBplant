using NetCDF
using Unitful
using Unitful: W, m, °C, hr, mol, K, s, J, Mg, kg, kPa, Pa

Threads.nthreads()

basepath = "/home/raf/Uni/Masters/NetCDF/"
load_netcdf(fpath, varname) = ncread(joinpath(basepath, fpath), varname)
load32(fpath, varname) = replace(load_netcdf(fpath, varname), -32768=>missing)
load64(fpath, varname) = replace(load_netcdf(fpath, varname), -2147483647=>missing)

# ncinfo(joinpath(basepath, fpath), "pot0cm_0pctShade/pot0cm_0pctShade_1999.nc") 


solr = load32("SOLR/SOLR_1999.nc", "SOLR") 
rh1cm = load32("RH1cm_0pctShade/RH1cm_0pctShade_1999.nc", "RH1cm")
rh120cm = load32("RH120cm/RH120cm_1999.nc", "RH120cm")
ta1cm = load32("Tair1cm_0pctShade/TA1cm_0pctShade_1999.nc", "TA1cm")
ta120cm = load32("Tair120cm/TA120cm_1999.nc", "TA120cm")
v1cm = load32("V1cm/V1cm_1999.nc", "V1cm")
v120cm = load32("V120cm/V120cm_1999.nc", "V120cm")
pot0cm = load64("pot0cm_0pctShade/pot0cm_0pctShade_1999.nc", "pot0cm")
pot2_5cm = load64("pot2.5cm_0pctShade/pot2.5cm_0pctShade_1999.nc", "pot2.5cm")
pot5cm = load64("pot5cm_0pctShade/pot5cm_0pctShade_1999.nc", "pot5cm")
pot10cm = load64("pot10cm_0pctShade/pot10cm_0pctShade_1999.nc", "pot10cm")
pot20cm = load64("pot20cm_0pctShade/pot20cm_0pctShade_1999.nc", "pot20cm")
pot30cm = load64("pot30cm_0pctShade/pot30cm_0pctShade_1999.nc", "pot30cm")
pot50cm = load64("pot50cm_0pctShade/pot50cm_0pctShade_1999.nc", "pot50cm")
pot100cm = load64("pot100cm_0pctShade/pot100cm_0pctShade_1999.nc", "pot100cm")
soil0cm = load32("soil5cm_0pctShade/soil5cm_0pctShade_1999.nc", "soil5cm")
soil2_5cm = load32("soil2.5cm_0pctShade/soil2.5cm_0pctShade_1999.nc", "soil2.5cm")
soil5cm = load32("soil0cm_0pctShade/soil0cm_0pctShade_1999.nc", "soil0cm")
soil10cm = load32("soil10cm_0pctShade/soil10cm_0pctShade_1999.nc", "soil10cm")
soil20cm = load32("soil20cm_0pctShade/soil20cm_0pctShade_1999.nc", "soil20cm")
soil30cm = load32("soil30cm_0pctShade/soil30cm_0pctShade_1999.nc", "soil30cm")
soil50cm = load32("soil50cm_0pctShade/soil50cm_0pctShade_1999.nc", "soil50cm")
soil100cm = load32("soil100cm_0pctShade/soil100cm_0pctShade_1999.nc", "soil100cm")

mask = solr[:,:,1]
tair = ta1cm, ta120cm
wind = v1cm, v120cm
relh = rh1cm, rh120cm
tsoil = soil0cm, soil2_5cm, soil5cm, soil10cm, soil20cm, soil30cm, soil50cm, soil100cm
pot = pot0cm, pot2_5cm, pot5cm, pot10cm, pot20cm, pot30cm, pot50cm, pot100cm
i = CartesianIndex(20, 30)

build_env(rad, tair, wind, relh, tsoil, pot, i::CartesianIndex) = begin
    Env(rad[i, :] .* 0.1 * W*m^-2, 
        map(x -> x[i, :] .* 0.1 .* °C .|> K, tair), 
        map(x -> x[i, :] .* 0.001, relh), 
        map(x -> x[i, :] .* 0.1 .* m*s^-1, wind), 
        map(x -> x[i, :] .* 0.1 * °C .|> K, tsoil),  
        map(x -> x[i, :] .* 0.1 * J * kg^-1 * 1Mg * m^-3 .|> kPa, pot)
       )  
end
