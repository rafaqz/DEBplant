using Unitful: hr

@load "long_env.jld"

plotly()

function solplot!(plt, model, u, envstart)
    model.dead[] = false
    model.environment_start[] = envstart
    tstop = 365 * 24 - 1
    prob = DiscreteProblem(model, u, (0, tstop))
    sol = solve(prob, FunctionMap(scale_by_time = true))
    shootvals = map(i -> sol[i][2], 1:tstop+1)
    rootvals = map(i -> sol[i][8] * -1, 1:tstop+1)
    rng = envstart.val:envstart.val+tstop
    plot!(plt, rng, shootvals, linecolor = :black)
    plot!(plt, rng, rootvals, linecolor = :grey)
end


model = uimodel
u = init_state(model)
envstart = 1000u"hr"
month_hours = 365.25 / 12 * 24hr 
global envstart_float = 1.0hr
plt = plot();
# for i in 1:27*12
for i in 1:100
    println("month: ", i)
    global envstart_float += month_hours
    envstart = round(Int, envstart_float.val) * hr
    solplot!(plt, model, u, envstart)
end
plot(plt, xlab="Time in hrs", ylab="Biomass in mols")
