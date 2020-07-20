using PlotThemes, LaTeXStrings, Dates

include(joinpath(dirname(@__FILE__), "load.jl"))

STARTMONTH = 7
MONTHS = [Dates.LOCALES["english"].months... Dates.LOCALES["english"].months...]
MONTH_HOURS = 365.25 / 12 * 24hr
YLIMS = (-3,4)
LIFESPAN = 4380
LINEWIDTH = 1.1
MULTIYEARLINEWIDTH = 0.2
DPI = 300
STARTYR = 2005
STOPYR = 2011
YEARS = STOPYR-STARTYR
SMALLFONT = 6
ALPHA = 1.0
ENVLABELS = (swp="Soil water\npotential (kPa)", vpd="Vapour press.\ndeficit", st="Soil\ntemperature (°C)")
ENVINCREMENTS = reshape([string.(Microclimate.LAYERINCREMENTS)...], 1, 8)
ENVRANGE = reshape([string.(Microclimate.LAYERRANGE)...], 1, 2)
SOILCOLORS = permutedims(reverse(get.(Ref(ColorSchemes.copper), 0.0:1/7:1)))

function plot_soiltemp(model, envname, envrange; kwargs...)
    st = soiltemperature(model.environment) .|> °C
    plot(ustrip.(st[envrange, 2:end]);
        color=SOILCOLORS,
        labels=ENVINCREMENTS[:, 2:end],
        xlab="",
        kwargs...
    )
end

function plot_swp(model, envname, envrange; kwargs...)
    swp = -1 .* ustrip(soilwaterpotential(model.environment))
    yticks = ([10^-1, 10^0, 10^1, 10^2, 10^3, 10^4], ["-10e-1", "-10e0", "-10e1", "-10e2", "-10e3", "-10e4"])
    plot(swp[envrange, 2:end];
        color=SOILCOLORS,
        labels=ENVINCREMENTS[:, 2:end],
        yscale=:log10,
        yflip=true,
        yticks=yticks,
        xlab="",
        kwargs...
    )
end

function plot_growth!(plt, model, u, envstart)
    model.dead[] = false
    model.environment_start[] = envstart
    tstop = LIFESPAN * 1.0hr
    prob = DiscreteProblem(model, u, (0.0hr, tstop))
    state = solve(prob, FunctionMap(scale_by_time = true))
    shootvals = map(i -> state[ustrip(round(typeof(1hr), i))][:VS] * 25/mol, 1hr:1hr:tstop+1hr)
    rootvals = map(i -> state[ustrip(round(typeof(1hr), i))][:VR] * -25/mol, 1hr:1hr:tstop+1hr)
    rng = envstart:1.0hr:envstart+tstop
    plot!(plt, rng, shootvals;
        linecolor=:black,
        label="Shoot",
        xlab="",
    )
    plot!(plt, rng, rootvals;
        linecolor=:grey,
        label="Root",
        xlab="",
    )
    return model, state
end

" Plot growth starting every moth in the available timespan"
function plot_months(model, u, envstart; kwargs...)
    growthplots = plot(;
        ylims=YLIMS,
        legend=:none,
        link=:x,
        kwargs...
    )
    for i in 1:YEARS*12
        println("month: ", i)
        envstart += MONTH_HOURS
        # Round to the start of a day
        envstart_hour = round(typeof(1hr), round(typeof(1d), envstart)) + 1hr
        plot_growth!(growthplots, model, u, envstart_hour)
    end
    growthplots
end

" Plot all starting states in the transect"
function plot_years(model, transect, u, envstart)
    locplots = []
    envrange = 1:1:size(radiation(model.environment), 1)
    xticks = (collect(0:8759*2:8760*YEARS), string.(collect(STARTYR:2:STOPYR)))

    # Plots for legends
    stateplot = plot(;
        grid=false,
        ylab="Plant structural mass",
        showaxis=false,
        xlims=(-2, -1),
        xlab="",
    )
    plot_growth!(stateplot, model, u, envstart)
    soilwaterplot = plot_swp(model, :t1, envrange;
         linewidth=MULTIYEARLINEWIDTH,
         xlims=(-100, -10),
         ylab="",
         grid=false,
         showaxis=false
    )
    soiltempplot = plot_soiltemp(model, :t1, envrange;
         linewidth=MULTIYEARLINEWIDTH,
         xlims=(-100, -10),
         ylab="",
         grid=false,
         showaxis=false
    )
    legends = plot(stateplot, soiltempplot, soilwaterplot;
        legend=:right,
        xticks=xticks,
        link=:x,
        xlab="",
        ylab="",
        margin=0px,
        bottom_margin=5px,
        left_margin=20px,
        layout=grid(3, 1),
    )
    for (i, envname) in enumerate(keys(transect))
        println(envname)
        model.environment = transect[envname]
        # Allometry Y intercept has to match seed size.
        println("using starting mass given in state for B0")
        if i == 1
            stateplot = plot_months(model, u, envstart; 
                ylab="Structural mass (g)",
                xshowaxis=false,
            )
            soiltempplot = plot_soiltemp(model, "", envrange; 
                ylab=ENVLABELS[:st],
                ylims=(-2, 50),
                xshowaxis=false,
                linewidth=MULTIYEARLINEWIDTH, 
                legend=false,
            )
            # Double plot this for flipped ylims range
            soilwaterplot = plot!(plot_swp(model, "", envrange; 
                ylab=ENVLABELS[:swp],
                xlab=envname,
                linewidth=MULTIYEARLINEWIDTH, 
                legend=false,
            ), (1, 20000.0))
            push!(locplots,
                  plot(stateplot, soiltempplot, soilwaterplot;
                      xticks=xticks,
                      link=:x,
                      margin=0px,
                      bottom_margin=5px,
                      left_margin=20px,
                      layout=grid(3, 1)
                  )
            )
        else
            stateplot = plot_months(model, u, envstart; 
                ylab="", 
                yshowaxis=false,
                xshowaxis=false,
            )
            soiltempplot = plot_soiltemp(model, envname, envrange; 
                ylab=" ",
                ylims=(-2, 50),
                yshowaxis=false, 
                xshowaxis=false,
                linewidth=MULTIYEARLINEWIDTH,
                legend=false,
            )
            soilwaterplot = plot!(plot_swp(model, envname, envrange; 
                ylab=" ",
                xlab=envname,
                yshowaxis=false, 
                linewidth=MULTIYEARLINEWIDTH,
                legend=false,
            ), (1, 20000.0))
            push!(locplots,
                plot(stateplot, soiltempplot, soilwaterplot;
                    xticks=xticks,
                    link=:x,
                    margin=0px,
                    bottom_margin=5px,
                    left_margin=-70px,
                    layout=grid(3, 1)
                )
            )
        end
    end
    plt = plot(locplots...,
        legends,
        size=(1200,700),
        dpi=DPI,
        layout=grid(1, length(locplots)+1, widths=envwidths(transect)),
    )
    plt
end

envwidths(envs) = begin
    n = length(keys(envs))
    [[0.9/n for i in 1:n]..., 0.1]
end

function calc_heights(microplots)
    num_microplots = length(microplots)
    h = 1/(num_microplots + 2)
    heights = vcat([2h], [h for x in 1:num_microplots])
    num_plots = 1 + num_microplots
    heights, num_plots
end

function plot_crossover(model, environment, u, envstart, months)
    tstop = round(typeof(1hr), months * MONTH_HOURS)
    plotstart = 1hr
    model.environment = environment
    model.environment_start[] = envstart
    model.dead[] = false
    tspan = ustrip.((plotstart, tstop))
    prob = DiscreteProblem(model, ustrip.(u), tspan)
    state = solve(prob, FunctionMap(scale_by_time = true))
    n = length(u) ÷ 2
    state = state'
    state[:, 4:6] = state[:, 4:6] * -1
    # Legends are separated out to keep them off the plots
    plts = plot_list(model, state, envstart, tstop, months, (1, ustrip(tstop)); 
        legend=false,
    )
    legends = plot_list(model, state, envstart, tstop, months, (-100, -10);
        legend=:right,
        xlims=(-100, -10),
        grid=false,
        showaxis=false,
        ylab=""
    )
    plot(plts, legends;
        layout=grid(1, 2; widths=[0.85, 0.15]),
        size=(1200,1500),
        dpi=DPI
    )
end

function plot_list(model, state, envstart, tstop, months, xlims; kwargs...)
    envrange = ustrip.(envstart:1hr:envstart+tstop)
    yticks = ([10^-1, 10^0, 10^1],
              string.(["-0.1", "-1", "-10"], Ref(" kPa")))
    xticks=(ustrip.(0hr:MONTH_HOURS:tstop), MONTHS[1+STARTMONTH:months+1+STARTMONTH])
    rad = radiation(model.environment)
    varrange = 1:size(state, 1)
    vars1 = model.records[1].vars
    vars2 = model.records[2].vars
    radplot = plot(rad[envrange];
        ylab="Radiation",
        xshowaxis=false,
        color=:black,
        kwargs...,
        legend=nothing,
    )
    stateplot = plot(state;
        linewidth=LINEWIDTH, color=[1 2 3],
        labels=reshape([STATELABELS...], 1, 6),
        ylab="State Variables\n(C/N mol)",
        xlab="",
        xshowaxis=false,
        kwargs...
    )
    assimplot = plot(ustrip.(hcat(model.records[1].J[2,1,varrange]));
        ylims=(0.0, 0.000017),
        ylabel="Assimilation\n(C-mol/hr)",
        xlabel="",
        xticks=xticks,
        xshowaxis=false,
        linewidth=LINEWIDTH,
        kwargs...,
        legend=false,
    )
    soiltempplot = plot_soiltemp(model, "", envrange; 
        ylab=ENVLABELS[:st],
        ylims=(5, 30),
        xshowaxis=false, 
        linewidth=LINEWIDTH, 
        kwargs...
    )
    soilwaterplot = plot_swp(model, "", envrange; 
        ylab=ENVLABELS[:swp],
        xshowaxis=false, 
        linewidth=LINEWIDTH, 
        kwargs...
    )
    tempcorrplot = plot(ustrip.(hcat(vars1.tempcorrection[varrange], vars2.tempcorrection[varrange]));
        ylabel="Temp. rate\ncorrection",
        xlabel="",
        linewidth=LINEWIDTH,
        labels=["Shoot" "Root"],
        xshowaxis=false,
        kwargs...
    )
    rateplot = plot(ustrip.(hcat(vars1.rate[varrange], vars2.rate[varrange]));
        ylabel="Growth rate\n(mol mol^-1 d^-1)",
        labels=["Shoot" "Root"],
        linewidth=LINEWIDTH,
        xshowaxis=false,
        kwargs...
    )
    depthplot = plot(ustrip.(hcat(vars1.height[varrange], vars2.height[varrange] .* -1));
        ylabel="Height/\nDepth (m)",
        labels=["Shoot" "Root"],
        linewidth=LINEWIDTH,
        xshowaxis=false,
        kwargs...
    )
    maxswpplot = plot(vars1.swp[varrange]; ylab="Available\nwater pot.",
        xticks=xticks,
        ylims=(-15, 0),
        linewidth=LINEWIDTH,
        labels=ENVRANGE,
        kwargs...,
        legend=false
    )
    plts = (stateplot, assimplot, rateplot, depthplot, soiltempplot, tempcorrplot, soilwaterplot, maxswpplot)
    x = 0.8/7
    plot(plts...;
        left_margin=30px,
        xlims=xlims,
        xticks=xticks,
        link=:x,
        layout=grid(length(plts), 1, heights=[0.2, x, x, x, x, x, x, x, x])
    )
end

function plot_byname(model, environment, u, envstart, months, plotids)
    tstop = round(typeof(1hr), months * MONTH_HOURS)
    envrange = ustrip.(envstart:1hr:envstart+tstop)
    vars1 = model.records[1].vars
    vars2 = model.records[2].vars
    plotstart = 1hr
    model.environment = environment
    model.environment_start[] = envstart
    model.dead[] = false
    tspan = ustrip.((plotstart, tstop))
    prob = DiscreteProblem(model, ustrip.(u), tspan)
    state = solve(prob, FunctionMap(scale_by_time = true))
    varrange = 1:size(state, 1)
    xticks=(ustrip.(0hr:MONTH_HOURS:tstop), MONTHS[1+STARTMONTH:months+1+STARTMONTH])
    yticks = (-0.05:0.05:0.2)
    st = soiltemperature(model.environment) .|> °C
    rh = relhumidity(model.environment)
    airtemp = airtemperature(model.environment)
    vpd = uconvert.(u"kPa", vapour_pressure_deficit.(airtemp, rh))
    n = length(u) ÷ 2
    mass = state' .* 25g
    mass[:, 4:6] = mass[:, 4:6] * -1
    massplot = plot(mass,
        linewidth=LINEWIDTH,
        color=[1 2 3],
        labels=reshape([STATELABELS...], 1, 6),
        ylabel="State",
        ytick=yticks,
        xticks=xticks,
        xlabel="",
        legend=:topleft,
    )
    assimplot = plot(ustrip.(hcat(model.records[1].J[2,1,:]));
        ylabel="Assimilation\n(C-mol/hr)",
        xlabel="",
        xticks=xticks,
        linewidth=LINEWIDTH,
        legend=nothing,
    )
    tempcorrplot = plot(ustrip.(hcat(vars1.tempcorrection[varrange], vars2.tempcorrection[varrange]));
        ylabel="Temp. rate\ncorrection",
        xlabel="",
        linewidth=LINEWIDTH,
        labels=["Shoot" "Root"],
        xshowaxis=false,
    )
    rateplot = plot(ustrip.(hcat(vars1.rate[varrange], vars2.rate[varrange]));
        ylabel="Growth rate\n(mol mol^-1 d^-1)",
        labels=["Shoot" "Root"],
        linewidth=LINEWIDTH,
        xshowaxis=false,
    )
    avswpplot = plot(hcat(model.records[1].vars.swp);
        ylabel="Available\nsoil water\npotential",
        xticks=xticks,
        linewidth=LINEWIDTH,
        labels=ENVRANGE,
        legend=nothing,
    )
    radplot = plot(radiation(environment); ylab="Radiation", color=:black, legend=:none, xlabel="Month")
    depthplot = plot(ustrip.(hcat(vars1.height[1:size(mass, 1)], vars2.height[1:size(mass, 1)] .* -1));
        ylabel="Height/\nDepth (m)",
        labels=["Shoot" "Root"],
        xticks=xticks,
        linewidth=LINEWIDTH,
        legend=:bottomleft,
    )
    soiltempplot = plot_soiltemp(model, "", envrange; 
        ylims=(-2, 50),
        xshowaxis=false, 
        linewidth=MULTIYEARLINEWIDTH, 
    )
    soilwaterplot = plot_swp(model, "", envrange; 
        xshowaxis=false, 
        linewidth=MULTIYEARLINEWIDTH
    )
    plots = (mass=massplot, assim=assimplot, avswp=avswpplot, depth=depthplot, 
             swp=soilwaterplot, st=soiltempplot, rad=radplot, 
             rate=rateplot, corr=tempcorrplot)
    selectedplots = [plots[id] for id in plotids]
    heights = assimheights(selectedplots, plotids)
    plot(selectedplots...;
        xlims=(1,ustrip(tstop)),
        left_margin=60px,
        link=:x,
        layout=grid(length(selectedplots), 1, heights=heights),
        size=(1300,700),
        dpi=DPI
    )
end

function assimheights(plts, plotids)
    l = length(plts)
    if plotids[1] == :mass
        h = 1/(l+1)
        [2h, [h for x in 1:l-1]...]
    else
        [1/l for x in plts]
    end
end
