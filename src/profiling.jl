# Model performance profiling

dir = "DEBSCRIPTS" in keys(ENV) ? ENV["DEBSCRIPTS"] : pwd()
include(joinpath(dir, "load.jl"))
environments, tspan = loadenvironments(dir)
environments[:controls] = MicroclimControl()

u = [1e-4, 1e-4, 1.0, 1e-4, 1e-4, 10.0]u"mol" # Initial value
du = fill(0.0u"mol/hr", 6); nothing
vars = (Vars(), Vars())
models = OrderedDict()
modeldir = joinpath(dir, "models")
include.(joinpath.(Ref(modeldir), readdir(modeldir)));
plant = models[:bb];

plant(du, u, nothing, 1u"hr");
prob = DiscreteProblem(plant, u, (0.0u"hr", 100.0u"hr"))
@time solve(prob, FunctionMap(scale_by_time = true));

function prof(plant, u; maxdepth=40, view=false)
    prob = DiscreteProblem(plant, u, (0.0u"hr", 100.0u"hr"))
    Profile.clear()
    solve(prob, FunctionMap(scale_by_time = true))
    @profile for i = 1:1000
        sol = solve(prob, FunctionMap(scale_by_time = true))
    end
    if view
        ProfileView.view()
    else
        Profile.print(maxdepth = maxdepth)
    end
end

using Profile, ProfileView
@profile show(Plant())
Profile.clear()
@profile for i = 1:1000000
    plant.dead[] = false
    plant(du, u, nothing, 1u"hr") 
end

#@profile MicroclimPoint(envgrid, CartesianIndex(56, 53))
ProfileView.view()

prof(plant, u; view=true)
