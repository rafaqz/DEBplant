using Revise
using DynamicEnergyBudgets, Unitful
using Flatten, OrdinaryDiffEq, Microclimate

du = [0.0 for i in 1:12]u"mol/hr"
u = [0.0, 1e-4, 0.0, 1e-4, 1e-4, 1e-4, 0.0, 1e-4, 0.0, 1e-4, 1e-4, 10.0]u"mol"
t = 0u"hr":1u"hr":8760u"hr"
du = [0.0 for i in 1:12]
u = [0.0, 1e-4, 0.0, 1e-4, 1e-4, 1e-4, 0.0, 1e-4, 0.0, 1e-4, 1e-4, 10.0]
t = 0:1:8760
# du = [0.0 for i in 1:18]u"mol/hr"
# u = [0.0, 1e-1, 0.0, 1e-1, 1e-1, 1e-1, 0.0, 1e-1, 0.0, 1e-1, 1e-1, 1e-1, 0.0, 1e-1, 0.0, 1e-1, 1e-1, 10.0]u"mol"
# u = [0.0, 1e-1, 0.0, 1e-1, 1e-1, 1e-1, 0.0, 1e-1, 0.0, 1e-1, 1e-1, 1.0]u"mol"
# const u = [0.0, 1e-1, 0.0, 1e-1, 1e-1, 1e-1, 0.0, 1e-1, 0.0, 1e-1, 1e-1, 1e-1, 0.0, 1e-1, 0.0, 1e-1, 1e-1, 1.0]u"mol"
env = environment

organism = models[:init]
organism = models[:maturity]
modelobs = Ref(DynamicEnergyBudgets.PlantCN(time=t, environment_start=1u"hr"));
organism = DynamicEnergyBudgets.PlantCN(time=t, environment_start=1u"hr");
organism = DynamicEnergyBudgets.PlantCN(environment=env, time=t, environment_start=1u"hr");
organism = DynamicEnergyBudgets.FvCBPlant(time=t);
organism = DynamicEnergyBudgets.FvCBPlant3(time=t);
organism = DynamicEnergyBudgets.FvCBPlant(environment=env, time=t);
organism(du, u, nothing, 10u"hr")
organism(du, u, nothing, 1)

length(organism.records[1].vars.rate)
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
