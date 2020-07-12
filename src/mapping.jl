include(joinpath(dirname(@__FILE__), "load.jl"))

using Statistics, Shapefile, GraphRecipes, StatsBase, Microclimate, NCDatasets, JLD2
using DynamicEnergyBudgets: dead

runsims(i, mask, model, u, envgrid, tspan, year) =
    if ismissing(mask) 
        missing
    else
        println(year, ": ", i)
        env = MicroclimPoint(envgrid, i)
        ismissing(env) && return missing
        model = @set model.environment = env
        # Round to the start of a day
        tstop = tspan - oneunit(tspan)
        u_sol = typeof(u)[]
        envstart = oneunit(MONTH_HOURS)
        # Run for each month
        while (envstart + tspan) < length(radiation(env)) * hr
            envstart_hour = round(typeof(1hr), round(typeof(1d), envstart))
            model.environment_start[] = envstart_hour
            # Run model in this environment
            model.dead[] = false
            model = set_allometry(model, u)
            sol = discrete_solve(model, u, tstop)
            push!(u_sol, dead(model) ? zero(sol) : sol)
            envstart += MONTH_HOURS
        end
        u_sol
    end

run_year(year, datapath, shade, model, skip) = begin
    years = year:year+1
    println("running $years")
    envgrid = load_grid(datapath, years, shade, skip)
    masklayer = airtemperature(envgrid)[1][1][:,:,1]
    tspan = 8759hr
    u = zeros(6)mol
    ulabelled = DimArray(u, X(Val(Tuple(Symbol.(STATELABELS)))))
    runsims.(CartesianIndices(masklayer), masklayer, Ref.((model, ulabelled, envgrid, tspan, year))...)
end

# envgrid = load_grid(datapath, years, shade, SKIPPED)
# masklayer = airtemperature(envgrid)[1][1][:,:,1]
# for i = CartesianIndices(masklayer)
    # println(i)
    # envs = MicroclimPoint(envgrid, i)
# end
#


combine_year(year) = begin 
    out = Array{Union{Missing,Float64},2}(undef, size(year)...)
    for i in CartesianIndices(year)
        out[i] = if ismissing(year[i]) 
            0.0 
        else
            trans(sum_structural_mass(year[i]))
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
sum_structural_mass(as) = maximum((ustrip.(A[:VS]) for A in as))

build_plot(data, name, legend) = begin
    data = rotl90(data) * 25
    hm = heatmap(long, lat, data; c=:tempo, title=name, clims=(0.0, 9.8), 
                 legend=legend, colorbar_title="Shoot structural mass (g)")
    plt = plot!(hm, shp.shapes; 
          xlim=(long[1]-1, long[end]+1), ylim=(lat[end]-1, lat[1]+1), 
          # xlab="Longitude", ylab="Latitude",
          color=:black, width=2, legend=false
         )
    longs = getindex.(Ref(long), [65, 60, 55])
    lats = getindex.(Ref(lat), [35, 35, 35])
    scatter!(plt, longs, lats; color=:white, markersize=2)
    annotate!(plt, longs, lats .+ 1, text.(["T1", "T2", "T3"], 7))
end

shape_plot(long, lat, points, labels, size, markersize) = begin
    plt = plot(shp.shapes; 
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
