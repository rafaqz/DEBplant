# include("NetCDF.jl")
model = uimodel;

using JLD2, LabelledArrays
using Unitful: mol

# function load_environment()
#     @load "model.jld"
#     model
# end

# @save "model.jld" model

solve_for_env(model, u, tstop) = begin
    prob = DiscreteProblem(model, u, (one(tstop), tstop) .* hr)
    solve(prob, FunctionMap(scale_by_time = true))
end

model.environment_start[] = oneunit(model.environment_start[])

statelabs = Symbol.((string.(STATE) .* "S"..., string.(STATE) .* "R"...))
initstate = [0.0, 1e-4, 0.0, 1e-4, 1e-4, 1e-4, 0.0, 1e-4, 0.0, 10, 2, 0.0] * mol
u = LVector{eltype(initstate),typeof(initstate),statelabs}(initstate)
tstop = length(solr[i, :])-1
sol = solve_for_env(model, u, tstop)
solutions = Union{Missing,typeof(sol.u[end])}[missing for i in CartesianIndices(mask)]

@time for i in CartesianIndices(mask)
    if !ismissing(mask[i]) 
        model.environment = build_env(solr, tair, wind, relh, tsoil, pot, i)
        model.dead[] = false
        sol = solve_for_env(model, u, tstop)
        solutions[i] = DynamicEnergyBudgets.dead(model) ? zero(sol.u[end]) : sol.u[end]
        # println("Dead?: ", model.dead[], sol.u[end])
    end
end

structure = map(s -> ismissing(s) ? s : s.VS, solutions) 
structure1 = replace(structure, missing => 0mol) ./ mol 
maximum(structure1)
heatmap(permutedims(structure1[1:end,end:-1:1]))
