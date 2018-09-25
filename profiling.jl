using Revise
using DynamicEnergyBudgets
using Unitful
using OrdinaryDiffEq

du = [0.0 for i in 1:12]u"mol/hr"
u = [0.0, 1e-4, 0.0, 1e-4, 1e-4, 1e-4, 0.0, 1e-4, 0.0, 1e-4, 1e-4, 10.0]u"mol"
t = 0u"hr":1u"hr":8760u"hr"
organism = DynamicEnergyBudgets.Plant(time=t);

using ProfileView
function prof(organism; maxdepth=40, view=false)
    prob = DiscreteProblem(organism, u, (0u"hr", 8750u"hr"))
    Profile.clear()
    solve(prob, FunctionMap(scale_by_time = true))
    for i = 1:100
        @profile sol = solve(prob, FunctionMap(scale_by_time = true))
    end
    if view
        ProfileView.view()
    else
        Profile.print(maxdepth = maxdepth)
    end
end

prof(organism; view = true)
# Profile.print(format=:flat)

