dir = "DEBSCRIPTS" in keys(ENV) ? ENV["DEBSCRIPTS"] : pwd()
include(joinpath(dir, "load.jl"))
include(joinpath(dir, "plantstates.jl"))
using Statistics

import Base: round
round(::Type{T}, x::Quantity) where {T<:Quantity} = T(round(typeof(one(T)), uconvert(unit(T), x).val))

using DynamicEnergyBudgets: dead

runsims(i, mask, model, smallseed, largeseed, plant, envgrid, tspan) =
    if ismissing(radiation(envgrid)[1][i,1]) 
        missing
    else
        println(i)
        env = MicroclimPoint(envgrid, i)
        model = @set model.environment = env
        # Round to the start of a day
        tstop = tspan - oneunit(tspan)
        smallseed_sol = typeof(smallseed)[]
        largeseed_sol = typeof(largeseed)[]
        plant_sol = typeof(plant)[]
        envstart = oneunit(MONTH_HOURS)
        # Run for each month
        while (envstart + tspan) < length(radiation(env)) * hr
            envstart_hour = round(typeof(1hr), round(typeof(1d), envstart))
            model.environment_start[] = envstart_hour
            # Run model in this environment
            model.dead[] = false
            model = set_allometry(model, smallseed)
            sol = discrete_solve(model, smallseed, tstop)
            push!(smallseed_sol, dead(model) ? zero(sol) : sol)
            model.dead[] = false
            model = set_allometry(model, largeseed)
            sol = discrete_solve(model, largeseed, tstop)
            push!(largeseed_sol, dead(model) ? zero(sol) : sol)
            model.dead[] = false
            sol = discrete_solve(model, plant, tstop)
            push!(plant_sol, dead(model) ? zero(sol) : sol)
            envstart += MONTH_HOURS
        end
        smallseed_sol, largeseed_sol, plant_sol 
    end

run_year(year, basepath, shade, model) = begin
    years = year:year+1
    println("running $years")
    envgrid = load_grid(basepath, years, shade, SKIPPED)
    masklayer = radiation(envgrid)[1][:,:,1]
    tspan = 8759hr
    u = zeros(12)mol
    ulabelled = LArray{LABELS}(u)
    runsims.(CartesianIndices(masklayer), masklayer, Ref.((model, smallseed, largeseed, plant, envgrid, tspan))...)
end

const SKIPPED = (:snowdepth, :soilwatercontent, :relhumidity, :windspeed) 
const LABELS = (:PS, :VS, :MS, :CS, :NS, :ES, :PR, :VR, :MR, :CR, :NR, :ER)
const MONTH_HOURS = 365.25 / 12 * 24hr

basepath = "/home/raf/Data/microclim"
years = 2001:2002
shade = 0
i = CartesianIndex(65,35)

# Import environments 
environments, tspan = loadenvironments(dir)

# Import models
models = OrderedDict()
modeldir = joinpath(dir, "models")
include.(joinpath.(Ref(modeldir), readdir(modeldir)));
model = deepcopy(models[:maturity]);
model.environment_start[] = oneunit(model.environment_start[])
yearly_outputs = run_year.(years, Ref(basepath), shade, Ref(model));

extract(yo, n) = [ismissing(yo[1][i]) ? 0.0 : sum_yo(yo, i, n) for i in ind]
sum_yo(yo, i, n) = mean((sum_VS(yo[x][i][n]) for x in eachindex(yo)))
sum_VS(a) = mean((ustrip(la.VS) for la in a))

ind = CartesianIndices(yearly_outputs[1])
smallseed_out = extract(yearly_outputs, 1)
largeseed_out = extract(yearly_outputs, 2) 
plant_out = extract(yearly_outputs, 3)

# Plot
gr()
heatmap(rotl90(plant_out), size=(600,400), dpi=100)
# savefig("plots/plant_aus.png")
heatmap(rotl90(smallseed_out), size=(600,400), dpi=100)
# savefig("plots/smallseed_aus.png")
heatmap(rotl90(largeseed_out), size=(600,400), dpi=100) 
# savefig("plots/largeseed_aus.png")

permutedims(plant_out, (2,1))
