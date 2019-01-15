using Revise
using NetCDF
using Unitful: W, m, °C, hr, mol, K
using DynamicEnergyBudgets: Env, PlantCN, STATE
using OrdinaryDiffEq
using LabelledArrays
using Plots

Threads.nthreads()

ncinfo("/home/raf/Uni/Masters/NetCDF/Tair1cm_0pctShade/TA1cm_0pctShade_1999.nc")
SOLR = replace(ncread("/home/raf/Uni/Masters/NetCDF/SOLR/SOLR_1999.nc", "SOLR"), -32768=>missing) 
TA1cm = replace(ncread("/home/raf/Uni/Masters/NetCDF/Tair1cm_0pctShade/TA1cm_0pctShade_1999.nc", "TA1cm"), -32768=>missing)
ncinfo("/home/raf/Uni/Masters/NetCDF/soil2.5cm_0pctShade/soil2.5cm_0pctShade_1999.nc")
soil0cm = replace(ncread("/home/raf/Uni/Masters/NetCDF/soil5cm_0pctShade/soil5cm_0pctShade_1999.nc", "soil5cm"), -32768=>missing)
soil2_5cm = replace(ncread("/home/raf/Uni/Masters/NetCDF/soil2.5cm_0pctShade/soil2.5cm_0pctShade_1999.nc", "soil2.5cm"), -32768=>missing)
soil5cm = replace(ncread("/home/raf/Uni/Masters/NetCDF/soil0cm_0pctShade/soil0cm_0pctShade_1999.nc", "soil0cm"), -32768=>missing)
soil10cm = replace(ncread("/home/raf/Uni/Masters/NetCDF/soil10cm_0pctShade/soil10cm_0pctShade_1999.nc", "soil10cm"), -32768=>missing)
soil20cm = replace(ncread("/home/raf/Uni/Masters/NetCDF/soil20cm_0pctShade/soil20cm_0pctShade_1999.nc", "soil20cm"), -32768=>missing)
soil30cm = replace(ncread("/home/raf/Uni/Masters/NetCDF/soil30cm_0pctShade/soil30cm_0pctShade_1999.nc", "soil30cm"), -32768=>missing)
soil50cm = replace(ncread("/home/raf/Uni/Masters/NetCDF/soil50cm_0pctShade/soil50cm_0pctShade_1999.nc", "soil50cm"), -32768=>missing)
soil100cm = replace(ncread("/home/raf/Uni/Masters/NetCDF/soil100cm_0pctShade/soil100cm_0pctShade_1999.nc", "soil100cm"), -32768=>missing)

mask = SOLR[:,:,1]
soiltemps = soil0cm, soil2_5cm, soil5cm, soil10cm, soil20cm, soil30cm, soil50cm, soil100cm
i = CartesianIndex(20, 20)

build_env(solr, ta, soil, i::CartesianIndex) = begin
    Env(solr[i, :] .* 0.1 * W*m^-2, 
        ta[i, :] .* 0.1 .* °C .|> K, 
        map(s -> s[i, :] .* 0.1 * °C .|> K, soil))  
end

env = build_env(SOLR, TA1cm, soiltemps, i)
E = typeof(env)
envdata = similar(mask, Union{Missing,E})

for i in CartesianIndices(mask)
    if !ismissing(mask[i]) 
        envdata[i] = build_env(SOLR, TA1cm, soiltemps, i)
    else
        envdata[i] = missing
    end
end

solve_for_env(env, u, tstop) = begin
    model = PlantCN(environment=env)
    prob = DiscreteProblem(model, u, (one(tstop), tstop) .* hr)
    solve(prob, FunctionMap(scale_by_time = true))
end

statelabels = Symbol.((string.(STATE) .* "S"..., string.(STATE) .* "R"...))
initstate = [0.0, 1e-4, 0.0, 1e-4, 1e-4, 1e-4, 0.0, 1e-4, 0.0, 10, 2, 0.0] * mol
u = LVector{eltype(initstate),typeof(initstate),statelabels}(initstate)
tstop = length(SOLR[i, :])
sol = solve_for_env(envdata[30, 30], u, tstop)
solutions = Union{Missing,typeof(sol.u[end])}[missingrfor i in CartesianIndices(mask)]

@time Threads.@threads for i in CartesianIndices(mask)
    if !ismissing(mask[i]) 
        @inbounds sol = solve_for_env(envdata[i], u, tstop)
        @inbounds solutions[i] = sol.u[end]
    end
end


structure = map(s -> ismissing(s) ? s : s.VS, solutions) 
heatmap(replace(structure, missing=>0))
heatmap(replace(structure, missing=>0mol) ./ mol)
