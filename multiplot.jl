dir = "DEBSCRIPTS" in keys(ENV) ? ENV["DEBSCRIPTS"] : pwd()
include(joinpath(dir, "load.jl"))
include(joinpath(dir, "plantstates.jl"))

const MONTH_HOURS = 365.25 / 12 * 24hr

import Base: round
round(::Type{T}, x::Quantity) where {T<:Quantity} = T(round(typeof(one(T)), uconvert(unit(T), x).val))

function solplot!(plt, model, u, envstart)
    model.dead[] = false
    model.environment_start[] = envstart
    tstop = 8759hr
    prob = DiscreteProblem(model, u, (0hr, tstop))
    sol = solve(prob, FunctionMap(scale_by_time = true))
    shootvals = map(i -> sol[ustrip(round(typeof(1hr), i))][2] * 25g/mol, 1hr:1hr:tstop+1hr)
    rootvals = map(i -> sol[ustrip(round(typeof(1hr), i))][8] * -25g/mol, 1hr:1hr:tstop+1hr)
    rng = envstart:1hr:envstart+tstop
    plot!(plt, rng, shootvals, linecolor = :black)
    plot!(plt, rng, rootvals, linecolor = :grey)
    return nothing
end

" Plot growth starting every moth in the available timespan"
function monthlyplot(title, model, u, envstart)
    println("Plotting: ", title)
    plt = plot()
    for i in 1:100
        println("month: ", i)
        envstart += month_hours
        # Round to the start of a day
        envstart_hour = round(typeof(1hr), round(typeof(1d), envstart)) + 1hr
        solplot!(plt, model, u, envstart_hour)
        yield()
    end
    plot(plt, plot_title=title, xlab="Time in hours", ylab="Structural biomass in grams (roots shown as negative)",
         legend=:none, size=(1200,800), dpi=100)
    savefig(title)
    return nothing
end

" Plot all starting states in all environments"
function monthlyplotall(model, environments, states, envstart)
    for envname in keys(environments)
        println(envname)
        for statename in keys(states)
            println(statename)
            model.environment = environments[envname]
            # Allometry Y intercept has to match seed size.
            model = set_allometry(model, states[statename])
            monthlyplot("plots/$(envname)_$(statename)", model, states[statename], envstart)
        end
    end
end

gr()

environments, _ = loadenvironments(dir)
models = OrderedDict()
modeldir = joinpath(dir, "models")
include.(joinpath.(Ref(modeldir), readdir(modeldir)));
envstart = 1.0hr

monthlyplotall(models[:maturity], environments, states, envstart)
