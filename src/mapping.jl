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
        max_state = typeof(u)[]
        envstart = oneunit(MONTH_HOURS)
        # Run for each month
        while (envstart + tspan) < length(radiation(env)) * hr
            envstart_hour = round(typeof(1hr), round(typeof(1d), envstart))
            model.environment_start[] = envstart_hour
            # Run model in this environment
            model.dead[] = false
            model = set_allometry(model, u)
            sol = max_state_solver(model, u, tstop)
            push!(max_state, sol)
            envstart += MONTH_HOURS
        end
        max_state
    end

# A simple manual solver is faster than DiffEq for this
function max_state_solver(model, u0, tstop)
    u = deepcopy(u0)
    max_ = deepcopy(u0)
    du = u ./ unit(tstop)
    for i = 2oneunit(tstop):oneunit(tstop):tstop
        model(du, u, nothing, i)
        model.dead[] && break
        u .+= du .* unit(tstop)
        max_ .= max.(max_, u)
    end
    # Return the maximum values
    max_
end

run_year(year, datapath, shade, model, skip) = begin
    years = year:year+1
    println("running $years")
    envgrid = load_grid(datapath, years, shade, skip)
    masklayer = airtemperature(envgrid)[1][1][:,:,1]
    tspan = 8759hr
    u = init_state(model)
    runsims.(CartesianIndices(masklayer), masklayer, Ref.((model, u, envgrid, tspan, year))...)
end


combine_year(year) = begin 
    out = Array{Union{Missing,Float64},2}(undef, size(year)...)
    for i in CartesianIndices(year)
        out[i] = if ismissing(year[i]) 
            missing
        else
            sum_structural_mass(year[i])
        end
    end
    out
end

combine_month(years) = begin 
    out = [zeros(Union{Missing,Float64}, size(years[1])) for x in 1:12]
    for year in years
        for i in CartesianIndices(year)
            for month in 1:12
                val = if ismissing(year[i]) || ismissing(year[i][month]) 
                    missing
                else
                    val = ustrip(year[i][month][:VS])
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
                missing
            else
                val = ustrip(year[i][2][month][:VS])
            end
            out[month][i] = max(out[month][i], val)
        end
    end
    out
end

sum_structural_mass(as) = maximum((ustrip.(A[:VS]) for A in as))

build_map(data, name, legend, lons, lats, t_lons, t_lats) = begin
    data = rotl90(data) 
    grams_per_mol = 25
    data .*= grams_per_mol
    hm = heatmap(lons, lats, data; 
        c=:tempo, 
        clims=(0.0, 9.8), 
        legend=legend, 
        title=name, 
        colorbar_title="Shoot structural mass (g)"
    )
    plt = plot!(hm, shp.shapes; 
        xlim=(lons[1]-1, lons[end]+1), 
        ylim=(lats[1]-1, lats[end]+1), 
        # xlab="Longitude", ylab="Latitude",
        color=:black, 
        width=2, 
        legend=false
    )
    scatter!(plt, t_lons, t_lats; 
        color=:red, 
        markerstrokecolor=:red, 
        markersize=3,
    )
    annotate!(plt, t_lons, t_lats .+ 1, text.(["T1", "T2", "T3"], 7))
end

shape_plot(lons, lats, t_lons, t_lats, labels, size, markersize) = begin
    plt = plot(shp.shapes; 
        xlim=(lons[1]-1, lons[end]+1), 
        ylim=(lats[1]-1, lats[end]+1), 
        # xlab="Longitude", ylab="Latitude",
        color=:black, 
        width=2, 
        legend=false, 
        size=size
    )
    scatter!(plt, t_lons, t_lats; color=:red, markersize=markersize)
    plot!(plt, t_lons, t_lats; color=:red, linestyle=:dot)
    annotate!(plt, t_lons, t_lats .+ 1, text.(labels, 7))
end
