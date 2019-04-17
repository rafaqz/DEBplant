dir = "DEBSCRIPTS" in keys(ENV) ? ENV["DEBSCRIPTS"] : pwd()
    rateplot = 
include(joinpath(dir, "load.jl"))
include(joinpath(dir, "plantstates.jl"))

const MONTH_HOURS = 365.25 / 12 * 24hr
YLIMS = (-450,650)
LIFESPAN = 4380

import Base: round
round(::Type{T}, x::Quantity) where {T<:Quantity} = T(round(typeof(one(T)), uconvert(unit(T), x).val))

function plot_microclim(model, envname, rnge = 1:1:size(radiation(model.environment), 1); 
                        ylabs=("Soil water\npotential", "Soil\ntemp.", "Radiation"), kwargs...)
    swp = -1 .* ustrip(soilwaterpotential(model.environment)) # .|> MPa
    st = soiltemperature(model.environment) .|> °C
    rad = radiation(model.environment)
    soilcolors = permutedims(reverse(get.(Ref(ColorSchemes.copper), 0.0:1/7:1)))
    # yticks = ["-10⁰ KPa", "-10¹ KPa", "-10² KPa", "-10³ KPa", "-10⁴ KPa"]
    # yticks = [-10^0u"kPa", -10^1u"kPa", -10^2u"kPa", -10^3u"kPa", -10^4u"kPa"]
    # yticks = [-10^0 -10^1 -10^2 -10^3 -10^4]
    soilwaterplot = plot(rnge * hr, swp[rnge, :]; yflip=true, scale=:log10, ylab=ylabs[1], 
                         color=soilcolors, legend=:none, xaxis=false, kwargs...)
    soiltempplot = plot(rnge * hr, st[rnge, :]; ylab=ylabs[2],
                        ylims=(-10, 80), color=soilcolors, legend=:none, xaxis=false, kwargs...)
    radplot = plot(rnge * hr, rad[rnge, :]; xlab=titlecase(string(envname)), ylab=ylabs[3], linealpha=0.8,
                   ylims=(0, 1100), legend=:none, color=:black, kwargs...)
    soilwaterplot, soiltempplot, radplot
end

function plot_sol!(plt, model, u, envstart)
    model.dead[] = false
    model.environment_start[] = envstart
    tstop = LIFESPAN * hr
    prob = DiscreteProblem(model, u, (0hr, tstop))
    sol = solve(prob, FunctionMap(scale_by_time = true))
    shootvals = map(i -> sol[ustrip(round(typeof(1hr), i))][2] * 25g/mol, 1hr:1hr:tstop+1hr)
    rootvals = map(i -> sol[ustrip(round(typeof(1hr), i))][8] * -25g/mol, 1hr:1hr:tstop+1hr)
    rng = envstart:1hr:envstart+tstop
    plot!(plt, rng, shootvals, linecolor=:black, xaxis=false)
    plot!(plt, rng, rootvals, linecolor=:grey, xaxis=false)
    return model, sol
end

" Plot growth starting every moth in the available timespan"
function plot_months(title, model, u, envstart; kwargs...)
    solplots = plot(; ylab=titlecase(string(title)), #ylims=YLIMS,
                    legend=:none, link=:x, kwargs...)

    for i in 1:119
        println("month: ", i)
        envstart += MONTH_HOURS
        # Round to the start of a day
        envstart_hour = round(typeof(1hr), round(typeof(1d), envstart)) + 1hr
        plot_sol!(solplots, model, u, envstart_hour)
    end
    solplots
end

" Plot all starting states in all environments"
function plot_years(model, environments, states, envstart)
    locplots = []
    for (i, envname) in enumerate(keys(environments))
        println(envname)
        stateplots = []

        for statename in keys(states)
            println("Plotting: ", statename)
            println(statename)
            model.environment = environments[envname]
            # Allometry Y intercept has to match seed size.
            model = set_allometry(model, states[statename])
            if i == 1
                push!(stateplots, plot_months(statename, model, states[statename],  envstart))
            else
                push!(stateplots, plot_months("", model, states[statename], envstart; yaxis=false))
            end
        end
        if i == 1
            microplots = plot_microclim(model, envname)
            push!(locplots, plot(stateplots..., microplots...;
                 xaxis=((0,8760*11), 0:30000:90000), link=:x, margin=0px, right_margin=0px,
                 layout=grid(6, 1, heights=[0.2666, 0.2666, 0.2666, 0.1, 0.05, 0.05])))
        else
            microplots = plot_microclim(model, envname; ylabs=("", "", ""), yaxis=false)
            push!(locplots, plot(stateplots..., microplots...;
                 xaxis=((0,8760*11), 0:30000:90000), link=:x, margin=0px, left_margin=-90px,
                 layout=grid(6, 1, heights=[0.2666, 0.2666, 0.2666, 0.1, 0.05, 0.05])))
        end
    end
    plt = plot(locplots..., size=(1200,1200), dpi=100, layout=grid(1, 3))
    savefig("plots/all")
    plt
end


function plot_single(model, environments, states, envstart)
    microplots = plot_microclim(model, "T1", 1:1:LIFESPAN; margin=0px)
    model.environment = environments[:t1]
    solplot = plot(yaxis=("Plant"), legend=:none, link=:x)
    model, sol = plot_sol!(solplot, model, states[:plant], envstart)
    rateplot = plot(hcat(model.records[1].vars.rate, model.records[2].vars.rate); yaxis=("Growth rate"))
    assimplot = plot(model.records[1].J[4,1,:]; yaxis=("Assimilation"))
    plt = plot(solplot, rateplot, assimplot, microplots...; xaxis=((0,LIFESPAN), 0:500:4000), link=:x,
         layout=grid(6, 1, heights=[0.3, 0.2, 0.15, 0.15, 0.1, 0.1]), size=(800,600), dpi=100)
    savefig("plots/single")
    plt
end


gr()
# pyplot()
# plotly()

environments, _ = loadenvironments(dir)
models = OrderedDict()
modeldir = joinpath(dir, "models")
include.(joinpath.(Ref(modeldir), readdir(modeldir)));
envstart = 1.0hr
model = models[:init]

plt = plot_single(model, environments, states, envstart);

allplot = plot_years(model, environments, states, envstart)

# plt = plot([1,2,3]hr, 1:3)
# plot!(plt, [1,2,3]hr, 4:6)
# plot(plt, plot(rand(6)), link=:x, layout=(2,1))
