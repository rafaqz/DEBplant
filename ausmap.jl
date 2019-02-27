dir = "DEBSCRIPTS" in keys(ENV) ? ENV["DEBSCRIPTS"] : pwd()
include(joinpath(dir, "load.jl"))
include(joinpath(dir, "plantstates.jl"))

using DynamicEnergyBudgets: dead

runsims(i, model, smallseed, largeseed, plant, envgrid, tstop) =
    if ismissing(masklayer[i]) 
        missing
    else
        println(i)
        model = @set model.environment = MicroclimPoint(envgrid, i)
        model.dead[] = false
        model = set_allometry(model, smallseed)
        sol = discrete_solve(model, smallseed, tstop)
        smallseed_sol = dead(model) ? zero(sol) : sol
        model.dead[] = false
        model = set_allometry(model, largeseed)
        sol = discrete_solve(model, largeseed, tstop)
        largeseed_sol = dead(model) ? zero(sol) : sol
        model.dead[] = false
        sol = discrete_solve(model, plant, tstop)
        plant_sol = dead(model) ? zero(sol) : sol
        # println(dead(model) ? "dead" : "alive")
        smallseed_sol, largeseed_sol, plant_sol 
    end

basepath = "/home/raf/Data/microclim"
years = 2001
shade = 0
i = CartesianIndex(65,35)

# Import environments 
environments, tspan = loadenvironments(dir)

# Import models
models = OrderedDict()
modeldir = joinpath(dir, "models")
include.(joinpath.(Ref(modeldir), readdir(modeldir)));
model = deepcopy(models[:maturity]);
model.environment_start[] = oneunit(model.environment_start[])

skipped = (:snowdepth, :soilwatercontent, :relhumidity, :windspeed) 
envgrid = load_grid(basepath, years, shade)
masklayer = radiation(envgrid)[1][:,:,1]

envpoint = MicroclimPoint(envgrid, i)
tstop = (length(radiation(envpoint))-1)hr
u = zeros(12)mol
labels = (:PS, :VS, :MS, :CS, :NS, :ES, :PR, :VR, :MR, :CR, :NR, :ER)
ulabelled = LArray{labels}(u)

env = MicroclimPoint(envgrid, i)
airtemperature(env)
model.records[1].vars.temp

sims = runsims.(CartesianIndices(masklayer), Ref.((model, smallseed, largeseed, plant, envgrid, tstop))...);

gr()

structure = replace(map(s -> ismissing(s) ? s : s[1][2], sims), missing => 0mol) ./ mol 
heatmap(permutedims(structure[1:end,end:-1:1]), size=(600,400), dpi=100)
savefig("plots/smallseed_aus.png")
structure = replace(map(s -> ismissing(s) ? s : s[2][2], sims), missing => 0mol) ./ mol 
heatmap(permutedims(structure[1:end,end:-1:1]), size=(600,400), dpi=100)
savefig("plots/largeseed_aus.png")
structure = replace(map(s -> ismissing(s) ? s : s[3][2], sims), missing => 0mol) ./ mol 
heatmap(permutedims(structure[1:end,end:-1:1]), size=(600,400), dpi=100)
savefig("plots/plant_aus.png")
