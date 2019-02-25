using Revise, LabelledArrays, DynamicEnergyBudgets, Photosynthesis, OrdinaryDiffEq, Unitful
using Plots, UnitfulPlots
using DynamicEnergyBudgets: allometry_pars, w_V
using Unitful: hr, d, mol, g, mg

dir = "DEBSCRIPTS" in keys(ENV) ? ENV["DEBSCRIPTS"] : pwd()
include(joinpath(dir, "load.jl"))

const month_hours = 365.25 / 12 * 24hr

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
        envstart = round(typeof(1hr), round(typeof(1d), envstart)) + 1hr
        solplot!(plt, model, u, envstart)
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

u = zeros(12)mol
ulabelled = LArray{STATEKEYS}(u)

smallseed = copy(ulabelled)
smallseed.VS = 1e-3mg / (25.0g/mol)
smallseed.CS = 1e-3mg / (25.0g/mol)
smallseed.NS = 1e-4mg / (25.0g/mol)
smallseed.VR = 1e-3mg / (25.0g/mol)
smallseed.CR = 1.0mg  / (25.0g/mol)
smallseed.NR = 0.05mg / (25.0g/mol)
smallseed
smallseed = [0.0, 1e-5, 0.0, 1e-5, 1e-5, 1e-5, 0.0, 1e-5, 0.0, 0.001, 0.00005, 0.0]mol

largeseed = copy(ulabelled)
largeseed.VS = 1e-1mg  / (25.0g/mol)
largeseed.CS = 1e-1mg  / (25.0g/mol)
largeseed.NS = 1e-1mg  / (25.0g/mol)
largeseed.VR = 1e-1mg  / (25.0g/mol)
largeseed.CR = 100.0mg / (25.0g/mol)
largeseed.NR = 5.0mg   / (25.0g/mol)
largeseed = [0.0, 1e-4, 0.0, 1e-4, 1e-4, 1e-4, 0.0, 1e-4, 0.0, 0.01, 0.0005, 0.0]mol

plant = copy(ulabelled)
plant.VS = 10g    / (25.0g/mol)
plant.CS = 10g    / (25.0g/mol)
plant.NS = 0.5g   / (25.0g/mol)
plant.VR = 5.0g   / (25.0g/mol)
plant.CR = 5.0mg  / (25.0g/mol)
plant.NR = 0.25mg / (25.0g/mol)
plant
plant = [0.0, 1e-2, 0.0, 1e-2, 1e-2, 1e-2, 0.0, 1e-2, 0.0, 1.0, 0.05, 0.0]mol

states = Dict(:smallseed => smallseed, :largeseed => largeseed, :plant => plant)

monthlyplotall(models[:maturity], environments, states, envstart)

module DocsCheck

using DocStringExtensions

@template TYPES =
    """
    $(TYPEDEF)
    $(DOCSTRING)
    """


macro documentbreaker(ex)
    :(Base.@__doc__ $(esc(ex)))
end

"Docs don't show typedef"
@documentbreaker struct BrokenDocs <: AbstractString 
    s::String
end


"Docs show typedef ^"
struct WorkingDocs <: AbstractString 
    s::String
end

end

DocsCheck.WorkingDocs
