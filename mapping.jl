include(joinpath(dir, "load.jl"))

using Statistics, Shapefile, GraphRecipes, StatsBase, Microclimate, NCDatasets, JLD2
using DynamicEnergyBudgets: dead

runsims(i, mask, model, largeseed, envgrid, tspan, year) =
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
        envstart = oneunit(MONTH_HOURS)
        # Run for each month
        while (envstart + tspan) < length(radiation(env)) * hr
            envstart_hour = round(typeof(1hr), round(typeof(1d), envstart))
            model.environment_start[] = envstart_hour
            # Run model in this environment
            model.dead[] = false
            model = set_allometry(model, largeseed)
            sol = discrete_solve(model, largeseed, tstop)
            push!(largeseed_sol, dead(model) ? zero(sol) : sol)
            envstart += MONTH_HOURS
        end
        largeseed_sol
    end

run_year(year, basepath, shade, model) = begin
    years = year:year+1
    println("running $years")
    envgrid = load_grid(basepath, years, shade, SKIPPED)
    masklayer = airtemperature(envgrid)[1][1][:,:,1]
    tspan = 8759hr
    u = zeros(12)mol
    ulabelled = LArray{LABELS}(u)
    runsims.(CartesianIndices(masklayer), masklayer, Ref.((model, largeseed, envgrid, tspan, year))...)
end

# envgrid = load_grid(basepath, years, shade, SKIPPED)
# masklayer = airtemperature(envgrid)[1][1][:,:,1]
# for i = CartesianIndices(masklayer)
    # println(i)
    # envs = MicroclimPoint(envgrid, i)
# end
#

#using NetCDF
const SKIPPED = (:snowdepth, :soilwatercontent) 
const LABELS = (:VS, :CS, :NS, :ES, :VR, :CR, :NR)
MONTH_HOURS = 365.25 / 12 * 24hr
basepath = "/home/raf/Data/microclim"
years = 2005:2010
shade = 0
i = CartesianIndex(65,35)
envgrid = load_grid(basepath, 2009:2009, shade, SKIPPED)
environments, _ = loadenvironments(dir)
# Import models
vars = (Vars(), Vars())
models = OrderedDict()
modeldir = joinpath(dir, "models")
include.(joinpath.(Ref(modeldir), readdir(modeldir)));
model = deepcopy(models[:bb]);
model.environment_start[] = oneunit(model.environment_start[])
scalingpath = joinpath(dir, "data/ausborder_polyline.shp")
shp = open(scalingpath) do io
    read(io, Shapefile.Handle)
end
# Get lat and long coordinates for plotting
radpath = joinpath(basepath, "SOLR/SOLR_2001.nc")
isfile(basepath)
long = NCDatasets.read(radpath, "longitude")
lat = NCDatasets.read(radpath, "latitude")
if isfile("data/yearly_outputs.jld")
    JLD2.@load "data/yearly_outputs.jld"
else
    # @time yearly_outputs = run_year.(years, Ref(basepath), shade, Ref(model));
    # JLD2.@save "data/yearly_outputs.jld" yearly_outputs 
end


combine_year(year) = begin 
    out = Array{Union{Missing,Float64},2}(undef, size(year)...)
    for i in CartesianIndices(year)
        out[i] = if ismissing(year[i]) 
            0.0 
        else
            trans(sum_VS(year[i]))
        end
    end
    out
end

combine_month(years) = begin 
    out = [zeros(Float64, size(years[1])...) for x in 1:12]
    for year in years
        for i in CartesianIndices(year)
            for month in 1:12
                val = if ismissing(year[i]) || ismissing(year[i][2][month]) 
                    0.0 
                else
                    val = trans(ustrip(year[i][month][:VS]))
                end
                out[month][i] = max(out[month][i], val)
            end
        end
    end
    out
end

extract_months(year) = begin
    out = [zeros(Float64, size(year)...) for x in 1:12]
    for i in CartesianIndices(year)
        for month in 1:12
            val = if ismissing(year[i]) || ismissing(year[i][2][month]) 
                0.0 
            else
                val = trans(ustrip(year[i][2][month][:VS]))
            end
            out[month][i] = max(out[month][i], val)
        end
    end
    out
end

trans(x) = x
sum_VS(a) = maximum((ustrip(la.VS) for la in a))

build_plot(data, name, legend) = begin
    data = rotl90(data) * 25
    hm = heatmap(long, lat, data; c=:tempo, title=name, clims=(0.0, 9.8), 
                 legend=legend, colorbar_title="Shoot structural mass (g)")
    plt = plot!(hm, shp.scalings; 
          xlim=(long[1]-1, long[end]+1), ylim=(lat[end]-1, lat[1]+1), 
          # xlab="Longitude", ylab="Latitude",
          color=:black, width=2, legend=false
         )
    longs = getindex.(Ref(long), [65, 60, 55])
    lats = getindex.(Ref(lat), [35, 35, 35])
    scatter!(plt, longs, lats; color=:white, markersize=2)
    annotate!(plt, longs, lats .+ 1, text.(["T1", "T2", "T3"], 7))
end

scaling_plot(long, lat, points, labels, size, markersize) = begin
    plt = plot(shp.scalings; 
          xlim=(long[1]-1, long[end]+1), ylim=(lat[end]-1, lat[1]+1), 
          # xlab="Longitude", ylab="Latitude",
          color=:black, width=2, legend=false, size=size
         )
    longs = points[1]
    lats = points[2]
    scatter!(plt, longs, lats; color=:red, markersize=markersize)
    plot!(plt, longs, lats; color=:red, linestyle=:dot)
    annotate!(plt, longs, lats .+ 1, text.(labels, 7))
end
