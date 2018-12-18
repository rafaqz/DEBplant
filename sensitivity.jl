using DiffEqSensitivity
using ForwardDiff

# Sensitivity
du = [0.0 for i in 1:12]
u = [0.0, 1e-2, 0.0, 1e-2, 1e-2, 1e-2, 0.0, 1e-2, 0.0, 1e-2, 1e-2, 10.0]
organism = DynamicEnergyBudgets.Plant();
values = flatten(Vector, organism)
fnames = fieldnameflatten(Vector, organism)

# vn = AxisArray(values, Axis{:parameters}(fnames))
# jacobianF = ForwardDiff.jacobian(p1->organism(du, u, p1, 1), flatten(Vector, organism))
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
