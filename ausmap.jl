dir = "DEBSCRIPTS" in keys(ENV) ? ENV["DEBSCRIPTS"] : pwd()
include(joinpath(dir, "load.jl"))
include(joinpath(dir, "plantstates.jl"))
using Statistics, Shapefile, GraphRecipes, StatsBase, Microclimate, NetCDF

import Base: round
round(::Type{T}, x::Quantity) where {T<:Quantity} = T(round(typeof(one(T)), uconvert(unit(T), x).val))

using DynamicEnergyBudgets: dead

runsims(i, mask, model, largeseed, plant, envgrid, tspan, year) =
    if ismissing(mask) 
        missing
    else
        println(year, ": ", i)
        env = MicroclimPoint(envgrid, i)
        ismissing(env) && return missing
        model = @set model.environment = env
        # Round to the start of a day
        tstop = tspan - oneunit(tspan)
        # smallseed_sol = typeof(smallseed)[]
        largeseed_sol = typeof(largeseed)[]
        plant_sol = typeof(plant)[]
        envstart = oneunit(MONTH_HOURS)
        # Run for each month
        while (envstart + tspan) < length(radiation(env)) * hr
            envstart_hour = round(typeof(1hr), round(typeof(1d), envstart))
            model.environment_start[] = envstart_hour
            # Run model in this environment
            model.dead[] = false
            model = set_allometry(model, largeseed)
            sol = discrete_solve(model, plant, tstop)
            push!(plant_sol, dead(model) ? zero(sol) : sol)
            model.dead[] = false
            # model = set_allometry(model, smallseed)
            # sol = discrete_solve(model, smallseed, tstop)
            # push!(smallseed_sol, dead(model) ? zero(sol) : sol)
            model.dead[] = false
            model = set_allometry(model, largeseed)
            sol = discrete_solve(model, largeseed, tstop)
            push!(largeseed_sol, dead(model) ? zero(sol) : sol)
            envstart += MONTH_HOURS
        end
        plant_sol, largeseed_sol
    end

run_year(year, basepath, shade, model) = begin
    years = year:year+1
    println("running $years")
    envgrid = load_grid(basepath, years, shade, SKIPPED)
    masklayer = airtemperature(envgrid)[1][1][:,:,1]
    tspan = 8759hr
    u = zeros(12)mol
    ulabelled = LArray{LABELS}(u)
    runsims.(CartesianIndices(masklayer), masklayer, Ref.((model, largeseed, plant, envgrid, tspan, year))...)
end

# envgrid = load_grid(basepath, years, shade, SKIPPED)
# masklayer = airtemperature(envgrid)[1][1][:,:,1]
# for i = CartesianIndices(masklayer)
    # println(i)
    # envs = MicroclimPoint(envgrid, i)
# end

const SKIPPED = (:snowdepth, :soilwatercontent) 
const LABELS = (:PS, :VS, :MS, :CS, :NS, :ES, :PR, :VR, :MR, :CR, :NR, :ER)
const MONTH_HOURS = 365.25 / 12 * 24hr
basepath = "/home/raf/Data/microclim"
years = 2001:2001
shade = 0
i = CartesianIndex(65,35)

environments, _ = loadenvironments(dir)

# Import models
models = OrderedDict()
modeldir = joinpath(dir, "models")
include.(joinpath.(Ref(modeldir), readdir(modeldir)));
model = deepcopy(models[:bbiso]);
model.environment_start[] = oneunit(model.environment_start[])
@time yearly_outputs = run_year.(years, Ref(basepath), shade, Ref(model));
shapepath = joinpath(basepath, "ausborder/ausborder_polyline.shp")

# Get lat and long coordinates for plotting
radpath = joinpath(basepath, "SOLR/SOLR_2001.nc")
long = ncread(radpath, "longitude")
lat = ncread(radpath, "latitude")


extract(yo, ind, n) = [ismissing(yo[1][i]) ? 0.0 : sum_yo(yo, i, n) for i in ind]
sum_yo(yo, i, n) = minimum((sum_VS(yo[x][i][n]) for x in eachindex(yo)))
sum_VS(a) = mean((ustrip(la.VS) for la in a))

# yearly_outputs = (yearly_outputs,)
# yearly_outputs = yearly_outputs[1]
ind = CartesianIndices(yearly_outputs[1])
sims = extract.(Ref(yearly_outputs), Ref(ind), (1,))


# Plot
gr()
# plotly()

shp = open(shapepath) do io
    read(io, Shapefile.Handle)
end

build_plot(data, name) = begin
    hm = heatmap(long, lat, rotl90(data); c=:tempo, legend=false, title=name)
    plot!(hm, shp.shapes; 
          xlim=(long[1], long[end]), ylim=(lat[end], lat[1]), 
          # xlab="Longitude", ylab="Latitude",
          color=:black, width=2, legend=false
         )
end

plts = build_plot.(sims, keys(states))
maps = plot(plts..., layout=(1,length(plts)), size=(1200,300), dpi=100)

# savefig("plots/spreadmaps.png")
