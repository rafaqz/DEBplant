using Revise
using JLD2, LabelledArrays
using Unitful: mol,m
using Distributed

@everywhere Microclimate, Unitful, DynamicEnergyBudgets

basepath = ENV["MICROCLIM"]
years = 2001
shade = 0
i = CartesianIndex(20,20)

@everywhere envgrid = load_grid(basepath, years, shade)
envpoint = MicroclimPoint(envgrid, CartesianIndex(20,20))



@everywhere model = uimodel;
model.environment_start[] = oneunit(model.environment_start[])
statelabs = Symbol.((string.(STATE) .* "S"..., string.(STATE) .* "R"...))
initstate = [0.0, 1e-4, 0.0, 1e-4, 1e-4, 1e-4, 0.0, 1e-4, 0.0, 10, 2, 0.0] * mol
u = LArray{statelabs}(initstate)

discrete_solve(model, u, tstop) = begin
    prob = DiscreteProblem(model, u, (one(tstop), tstop) .* hr)
    solve(prob, FunctionMap(scale_by_time = true))
end


masklayer = snowdepth(envgrid)[1][:,:,1]
tstop = length(radiation(envpoint))-1
sol = discrete_solve(model, u, tstop)

@everywhere begin
    tstop = length(radiation(envpoint))-1
    u = zeros(12)mol
    labels = (:PS, :VS, :MS, :CS, :NS, :ES, :PR, :VR, :MR, :CR, :NR, :ER)
    ulabelled = LArray{labels}(u)

    seed = copy(ulabelled)
    seed.VS = 1e-4mol 
    seed.CS = 1e-4mol 
    seed.NS = 1e-4mol 
    seed.VR = 1e-4mol 
    seed.CR = 0.01mol
    seed.NR = 0.0005mol

    plant = copy(ulabelled)
    plant.VS = 10mol 
    plant.CS = 5mol
    plant.NS = 0.5mol
    plant.VR = 10mol
    plant.CR = 5mol
    plant.NR = 0.5mol
end

u = seed

envstart = 1.0hr

function runsims(i, model, plant, seed, ulabelled, envgrid) 
    if !ismissing(masklayer[i]) 
        println(i)
        model.environment = MicroclimPoint(envgrid, i)
        model.dead[] = false
        try
            sol = discrete_solve(model, seed, tstop)
        catch
            seed_solutions = copy(ulabelled)
        finally
            seed_solutions = DynamicEnergyBudgets.dead(model) ? zero(sol.u[end]) : sol.u[end]
        end
        try
            sol = discrete_solve(model, plant, tstop)
        catch
            plant_solutions = copy(ulabelled)
        finally
            plant_solutions = DynamicEnergyBudgets.dead(model) ? zero(sol.u[end]) : sol.u[end]
        end
        # println("Dead?: ", model.dead[], sol.u[end])
    end
    seed_solutions, plant_solutions 
end

@distributed for i in CartesianIndices(masklayer) 
    runsims(i, model, plant, seed, ulabelled, envgrid) 
end

structure = map(s -> ismissing(s) ? s : s.VS, solutions) 
structure1 = replace(structure, missing => 0mol) ./ mol 
maximum(structure1)
heatmap(permutedims(structure1[1:end,end:-1:1]))
