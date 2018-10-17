using Revise
using DynamicEnergyBudgets
using Unitful
using Flatten
using OrdinaryDiffEq
using DifferentialEquations

using Plots

u0 = [10.0, 10.0, 10.0, 10.0]
prob = ODEProblem(f, u0, (0.0, 100.0))
sol = solve(prob)

# using Photosynthesis
# using TypedTables
# using AxisArrays
# using DiffEqSensitivity
# using ForwardDiff
# using JuliaDB
using IterableTables, DataFrames
using UnitfulPlots
using InteractBulma, Blink, TableWidgets, Plots, InteractBase, WebIO, Observables, Plots, CSSUtil

include("../DynamicEnergyBudgets/scratch/hide.jl")

du = [0.0 for i in 1:12]
u = [0.0, 1e-2, 0.0, 1e-2, 1e-2, 1e-2, 0.0, 1e-2, 0.0, 1e-2, 1e-2, 10.0]
organism = DynamicEnergyBudgets.Plant();
values = flatten(Vector, organism)
fnames = metaflatten(Vector, organism, fieldname_meta)
vn = AxisArray(values, Axis{:parameters}(fnames))
jacobianF = ForwardDiff.jacobian(p1->organism(du, u, p1, 1), flatten(Vector, organism))

organism = DynamicEnergyBudgets.Plant();
# organism = DynamicEnergyBudgets.FvCBPlant();
params = organism
def = Flatten.flatten(Vector, params)
ran = Flatten.metaflatten(Vector, params, DynamicEnergyBudgets.range)
pri = metaflatten(Vector, params, DynamicEnergyBudgets.prior)
fnames = metaflatten(Vector, params, fieldname_meta)
unit = metaflatten(Vector, params, DynamicEnergyBudgets.units)
form = metaflatten(Vector, params, fieldparent_meta)
lab = metaflatten(Vector, params, DynamicEnergyBudgets.label)
t = DataFrame(Formulation=form, Name=fnames, Default=def, Units=unit, Prior=pri, Label=lab)

# collect(zip(reshape(sum(jacobianF,1),length(p)), flatten_fieldnames(organism)))
# collect(zip(p, flatten_fieldnames(organism)))

prob = ODELocalSensitivityProblem(organism,u,(0, 50), vn)
sol = solve(prob, FunctionMap(scale_by_time = true))
x, dp = extract_local_sensitivities(sol)
x
ndp = AxisArray(dp, Axis{:parameters}(names))
da = ndp[:k_EN]
save("../DynamicEnergyBudgets/scratch/sensitivity.jld", "sensitivity", sol)


n = [flatten_fieldnames(organism)...]
da = ndp[:j_E_mai]
plot(da')
STATE = DynamicEnergyBudgets.STATE
labels = reshape(string.(vcat(STATE,STATE,STATE)), 1, 18)


du = [0.0 for i in 1:12]u"mol/hr"
u = [0.0, 1e-4, 0.0, 1e-4, 1e-4, 1e-4, 0.0, 1e-4, 0.0, 1e-4, 1e-4, 10.0]u"mol"
t = 0u"hr":1u"hr":8760u"hr"
# du = [0.0 for i in 1:18]u"mol/hr"
# u = [0.0, 1e-1, 0.0, 1e-1, 1e-1, 1e-1, 0.0, 1e-1, 0.0, 1e-1, 1e-1, 1e-1, 0.0, 1e-1, 0.0, 1e-1, 1e-1, 10.0]u"mol"
# u = [0.0, 1e-1, 0.0, 1e-1, 1e-1, 1e-1, 0.0, 1e-1, 0.0, 1e-1, 1e-1, 1.0]u"mol"
# const u = [0.0, 1e-1, 0.0, 1e-1, 1e-1, 1e-1, 0.0, 1e-1, 0.0, 1e-1, 1e-1, 1e-1, 0.0, 1e-1, 0.0, 1e-1, 1e-1, 1.0]u"mol"

organism = DynamicEnergyBudgets.Plant(time=t);
typeof(reconstruct(organism, flatten(Tuple, organism)))
typeof(reconstruct(organism, flatten(Vector, organism)))
organism = DynamicEnergyBudgets.Plant(environment=env2, time=t);
organism = DynamicEnergyBudgets.FvCBPlant(time=t);
organism = DynamicEnergyBudgets.FvCBPlant3(time=t);
organism = DynamicEnergyBudgets.FvCBPlant(environment=env2, time=t);
organism = DynamicEnergyBudgets.Plant(time=t);
organism = DynamicEnergyBudgets.PlantCN(time=t);
organism(du, u, nothing, 1u"hr")



prob = DiscreteProblem(organism, u, (0u"hr", 6000u"hr"))
sol = solve(prob, FunctionMap(scale_by_time = true))

plot(sol)
p = plot(sol.t, sol', size=(1000, 700))

# solve(prob, FunctionMap(scale_by_time = true))
using BenchmarkTools
Profile.clear()
@time for i = 1:100 solve(prob, FunctionMap(scale_by_time = true)) end
@btime $organism($du, $u, $vn, 1u"hr");
@btime $organism($du, $u, nothing, 1u"hr");


using JLD
using Microclimate
using IterableTables
using DataFrames
using TypedTables
# using IndexedTables
environment = load("../DynamicEnergyBudgets/scratch/environment.jld")["environment"]
# environment = nichemap_global("Adelaide", years=10)

env2 = Microclimate.MicroclimateTable(
  Table(environment.soil),
  Table(environment.shadsoil),
  Table(environment.metout),
  Table(environment.shadmet),
  Table(environment.soilmoist),
  Table(environment.shadmoist),
  Table(environment.humid),
  Table(environment.shadhumid),
  Table(environment.soilpot),
  Table(environment.shadpot),
  Table(environment.plant),
  Table(environment.shadplant),
  environment.RAINFALL,
  environment.dim,
  environment.ALTT,
  environment.REFL,
  environment.MAXSHADES,
  environment.longlat,
  environment.nyears,
  environment.timeinterval,
  environment.minshade,
  environment.maxshade,
  environment.DEP,
)
