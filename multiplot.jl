using JLD2, LabelledArrays, DynamicEnergyBudgets
using UnitfulPlots
using Unitful: hr, d

const month_hours = 365.25 / 12 * 24hr 

import Base: round
round(::Type{T}, x::Quantity) where {T<:Quantity} = T(($f)(typeof(one(T)), uconvert(unit(T), x).val))

function solplot!(plt, model, u, envstart)
    model.dead[] = false
    model.environment_start[] = envstart
    tstop = 365 * 24hr
    prob = DiscreteProblem(model, u, (0hr, tstop))
    local sol
    try
        sol = solve(prob, FunctionMap(scale_by_time = true))
    catch
        return nothing
    end
    shootvals = map(i -> sol[ustrip(round(typeof(1hr), i))][2] * 25g/mol, 1hr:1hr:tstop+1hr)
    rootvals = map(i -> sol[ustrip(round(typeof(1hr), i))][8] * -25g/mol, 1hr:1hr:tstop+1hr)
    rng = envstart:1hr:envstart+tstop
    plot!(plt, rng, shootvals, linecolor = :black)
    plot!(plt, rng, rootvals, linecolor = :grey)
    return nothing
end

function multiplot(title, model, u, envstart)
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
         legend=:none, size=(1200,800), dpi=300)
    savefig(title)
    return nothing
end

dir = "DEBSCRIPTS" in keys(ENV) ? ENV["DEBSCRIPTS"] : pwd()
include(joinpath(dir, "app.jl"))

environments, _ = loadenvironments(dir)
models = loadsavedmodels(dir)
model = models[:init]

u = zeros(12)g
labels = (:PS, :VS, :MS, :CS, :NS, :ES, :PR, :VR, :MR, :CR, :NR, :ER)
ulabelled = LVector{eltype(u),typeof(u),labels}(u)

small_seed = copy(ulabelled)
small_seed.VS = 1e-3mg * 25.0g/mol 
small_seed.CS = 1e-3mg * 25.0g/mol 
small_seed.NS = 1e-3mg * 25.0g/mol 
small_seed.VR = 1e-3mg * 25.0g/mol 
small_seed.CR = 1.0mg * 25.0g/mol
small_seed.NR = 0.05mg * 25.0g/mol
small_seed

large_seed = copy(ulabelled)
large_seed.VS = 1e-1mg * 25.0g/mol 
large_seed.CS = 1e-1mg * 25.0g/mol 
large_seed.NS = 1e-1mg * 25.0g/mol 
large_seed.VR = 1e-1mg * 25.0g/mol 
large_seed.CR = 100.0mg * 25.0g/mol
large_seed.NR = 5.0mg * 25.0g/mol
large_seed

plant = copy(ulabelled)
plant.VS = 10g * 25.0g/mol 
plant.CS = 10g * 25.0g/mol 
plant.NS = 0.5g * 25.0g/mol 
plant.VR = 5.0g * 25.0g/mol 
plant.CR = 5.0mg * 25.0g/mol
plant.NR = 0.25mg * 25.0g/mol
plant

u = seed

envstart = 1.0hr

model.environment = tas
multiplot("Tasmania_seed", model, seed, envstart)
multiplot("Tasmania_plant", model, plant, envstart)
model.environment = desert
multiplot("Desert_seed", model, seed, envstart)
multiplot("Desert_plant", model, plant, envstart)
model.environment = qld
multiplot("QLD_seed", model, seed, envstart)
multiplot("QLD_plant", model, plant, envstart)
