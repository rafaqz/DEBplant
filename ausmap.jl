using Revise
using JLD2, LabelledArrays, Microclimate
using Unitful: mol,m

basepath = ENV["MICROCLIM"]
years = 2001:2002
shade = 0
i = CartesianIndex(20,20)

envgrid = load_grid(basepath, years, shade)
envpoint = MicroclimPoint(envgrid, CartesianIndex(20,20))


model.environment_start[] = oneunit(model.environment_start[])
statelabs = Symbol.((string.(STATE) .* "S"..., string.(STATE) .* "R"...))
initstate = [0.0, 1e-4, 0.0, 1e-4, 1e-4, 1e-4, 0.0, 1e-4, 0.0, 10, 2, 0.0] * mol
u = LVector{eltype(initstate),typeof(initstate),statelabs}(initstate)

model = uimodel;

discrete_solve(model, u, tstop) = begin
    prob = DiscreteProblem(model, u, (one(tstop), tstop) .* hr)
    solve(prob, FunctionMap(scale_by_time = true))
end


tstop = length(solr[i, :])-1
sol = discrete_solve(model, u, tstop)
seed_solutions = Union{Missing,typeof(sol.u[end])}[missing for i in CartesianIndices(mask)]
solutions = Union{Missing,typeof(sol.u[end])}[missing for i in CartesianIndices(mask)]

u = zeros(12)mol
labels = (:PS, :VS, :MS, :CS, :NS, :ES, :PR, :VR, :MR, :CR, :NR, :ER)
ulabelled = LVector{eltype(u),typeof(u),labels}(u)

seed = copy(ulabelled)
seed.VS = 1e-4mol 
seed.CS = 1e-4mol 
seed.NS = 1e-4mol 
seed.VR = 1e-4mol 
seed.CR = 0.01mol
seed.NR = 0.0005mol
seed

plant = copy(ulabelled)
plant.VS = 10mol 
plant.CS = 5mol
plant.NS = 0.5mol
plant.VR = 10mol
plant.CR = 5mol
plant.NR = 0.5mol
plant

u = seed

envstart = 1.0hr

@time for i in CartesianIndices(mask)
    if !ismissing(mask[i]) 
        model.environment = MicroclimPoint(envgrid, i)
        model.dead[] = false
        sol = discrete_solve(model, seed, tstop)
        seed_solutions[i] = DynamicEnergyBudgets.dead(model) ? zero(sol.u[end]) : sol.u[end]
        sol = discrete_solve(model, plant, tstop)
        plant_solutions[i] = DynamicEnergyBudgets.dead(model) ? zero(sol.u[end]) : sol.u[end]
        # println("Dead?: ", model.dead[], sol.u[end])
    end
end

structure = map(s -> ismissing(s) ? s : s.VS, solutions) 
structure1 = replace(structure, missing => 0mol) ./ mol 
maximum(structure1)
heatmap(permutedims(structure1[1:end,end:-1:1]))
