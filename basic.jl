dir = "DEBSCRIPTS" in keys(ENV) ? ENV["DEBSCRIPTS"] : pwd()
include(joinpath(dir, "util/hide.jl"))
include(joinpath(dir, "load.jl"))

du = [0.0 for i in 1:6]u"mol/hr"
u = [0.0, 1e-4, 0.0, 1e-4, 1e-4, 1e-4, 0.0, 1e-4, 0.0, 1e-4, 1e-4, 10.0]u"mol"
u = [1e-4, 1e-4, 1e-4, 1e-4, 1e-4, 10.0]u"mol"
t = 0u"hr":1u"hr":8760u"hr"
# du = [0.0 for i in 1:12]
# u = [0.0, 1e-4, 0.0, 1e-4, 1e-4, 1e-4, 0.0, 1e-4, 0.0, 1e-4, 1e-4, 10.0]
# t = 0:1:8760

# Import environments 
environments, tspan = loadenvironments(dir)
env = first(values(environments))

# Import models
models = OrderedDict()
modeldir = joinpath(dir, "models")
include.(joinpath.(Ref(modeldir), readdir(modeldir)));

organism = models[:init];
organism = models[:maturity];
organism = models[:bbiso];
organism = app.savedmodel
v = organism.records[1].vars
organism.environment = environments[:t1];
organism = DynamicEnergyBudgets.PlantCN(time=t, environment_start=1u"hr");
organism = DynamicEnergyBudgets.PlantCN(environment=env, time=t, environment_start=1u"hr");
organism = DynamicEnergyBudgets.FvCBPlant(time=t);
organism = DynamicEnergyBudgets.FvCBPlant(environment=env, time=t);
organism(du, u, nothing, 10u"hr")
organism(du, u, nothing, 1)
fieldnames(typeof(organism.params[1].assimilation_pars))

organism

organism.records[1].vars.height

organs = define_organs(organism, 1hr)
o = organs[1]
ux = DynamicEnergyBudgets.split_state(organs, u)
DynamicEnergyBudgets.mass(o, ux[1])
DynamicEnergyBudgets.update_height!(DynamicEnergyBudgets.Allometry(), o, ux[1])
DynamicEnergyBudgets.update_height!(DynamicEnergyBudgets.FixedAllometry(), o, ux[1])
DynamicEnergyBudgets.update_height!(DynamicEnergyBudgets.SqrtAllometry(), o, ux[1])

prob = DiscreteProblem(organism, u, (0u"hr", 1000u"hr"))
sol = solve(prob, FunctionMap(scale_by_time = true))

plot(sol)
s = sol'[:, 7:12] .*= -1
plot(sol.t, sol', size=(1000, 700))

using BenchmarkTools
Profile.clear()
@time for i = 1:100 solve(prob, FunctionMap(scale_by_time = true)) end
@btime $organism($du, $u, $vn, 1u"hr");
@btime $organism($du, $u, nothing, 1u"hr");

def = Flatten.flatten(Vector, organism)
ran = Flatten.metaflatten(Vector, organism, DynamicEnergyBudgets.range)
pri = metaflatten(Vector, organism, DynamicEnergyBudgets.prior)
fnames = metaflatten(Vector, organism, fieldname_meta)
unit = metaflatten(Vector, organism, DynamicEnergyBudgets.units)
form = metaflatten(Vector, organism, fieldparent_meta)
lab = metaflatten(Vector, organism, DynamicEnergyBudgets.label)
t = DataFrame(Formulation=form, Name=fnames, Default=def, Units=unit, Prior=pri, Label=lab)
