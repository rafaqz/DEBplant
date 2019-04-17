dir = "DEBSCRIPTS" in keys(ENV) ? ENV["DEBSCRIPTS"] : pwd()
include(joinpath(dir, "load.jl"))
include(joinpath(dir, "plantstates.jl"))

function plotsol!(plt, model, u, envstart)
    model.dead[] = false
    model.environment_start[] = envstart
    tstop = 8759hr
    prob = DiscreteProblem(model, u, (0hr, tstop))
    sol = solve(prob, FunctionMap(scale_by_time = true))
    shootvals = map(i -> sol[ustrip(round(typeof(1hr), i))][2] * 25g/mol, 1hr:1hr:tstop+1hr)
    rootvals = map(i -> sol[ustrip(round(typeof(1hr), i))][8] * -25g/mol, 1hr:1hr:tstop+1hr)
    rng = envstart:1hr:envstart+tstop
    plot!(plt, rng, shootvals, linecolor=:black)
    plot!(plt, rng, rootvals, linecolor=:grey)
    return nothing
end
