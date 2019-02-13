using JLD2, LabelledArrays, DynamicEnergyBudgets
using UnitfulPlots
using Unitful: hr, d

const month_hours = 365.25 / 12 * 24hr 

import Base: round
# @eval ($f)(::Type{T}, x::Quantity) where {T<:Quantity} = T(($f)(typeof(one(T)), uconvert(unit(T), x).val))

@load "locations.jld"

plotly()

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

gr()

model = uimodel;

u = zeros(12)mol
labels = (:PS, :VS, :MS, :CS, :NS, :ES, :PR, :VR, :MR, :CR, :NR, :ER)
ulabelled = LVector{eltype(u),typeof(u),labels}(u)

seed = copy(ulabelled)
seed.VS = 1e-4mol 
seed.CS = 1e-4mol 
seed.NS = 1e-4mol 
seed.VR = 1e-4mol 
seed.CR = 0.01mol
seed.NR = 0.0005mol
seed

plant = copy(ulabelled)
plant.VS = 10mol 
plant.CS = 5mol
plant.NS = 0.5mol
plant.VR = 10mol
plant.CR = 5mol
plant.NR = 0.5mol
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
