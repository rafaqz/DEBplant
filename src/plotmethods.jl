using PlotThemes, LaTeXStrings, Dates

include(joinpath(dirname(@__FILE__), "load.jl"))

STARTMONTH = 7
MONTHS = [Dates.LOCALES["english"].months... Dates.LOCALES["english"].months...]
MONTH_HOURS = 365.25 / 12 * 24hr
YLIMS = (-3,4)
LIFESPAN = 4380
LINEWIDTH = 1.3
MULTIYEARLINEWIDTH = 0.2
DPI = 150
STARTYR = 2005
STOPYR = 2011
YEARS = STOPYR-STARTYR
SMALLFONT = 6
ALPHA = 1.0
ENVLABELS = (swp="Soil water\npotential", vpd="Vapour press.\ndeficit", st="Soil\ntemperature")
ENVINCREMENTS = reshape([string.(Microclimate.LAYERINCREMENTS)...], 1, 8)
ENVRANGE = reshape([string.(Microclimate.LAYERRANGE)...], 1, 2)

function plot_microclim(model, envname, rnge = 1:1:size(radiation(model.environment), 1);
                        ylabs=ENVLABELS, kwargs...)
    swp = -1 .* ustrip(soilwaterpotential(model.environment)) # .|> MPa
    st = soiltemperature(model.environment) .|> °C
    rh = relhumidity(model.environment)
    airtemp = airtemperature(model.environment)
    vpd = uconvert.(u"kPa", vapour_pressure_deficit.(airtemp, rh))
    soilcolors = permutedims(reverse(get.(Ref(ColorSchemes.copper), 0.0:1/7:1)))
    yticks = ([10^-1, 10^0, 10^1, 10^2, 10^3, 10^4],
              string.(["-10e-1", "-10e0", "-10e1", "-10e2", "-10e3", "-10e4"], Ref(" kPa")))
    soilwaterplot = plot(rnge * hr, swp[rnge, 2:end];
        yscale=:log10,
        ylab=ylabs[:swp],
        yflip=true,
        yticks=yticks,
        color=soilcolors,
        labels=ENVINCREMENTS[:, 2:end],
        kwargs...
    )
    # Make the plot full range keeping log scale. But not for the legend...
    if !(:showaxis in keys(kwargs))
        plot!(soilwaterplot, (1, 20000.0))
        plot!(soilwaterplot, (1, 0.1))
    end
    soiltempplot = plot(st[rnge, 2:end];
        ylab=ylabs[:st],
        ylims=(-2, 50),
        color=soilcolors,
        xlab=envname,
        labels=ENVINCREMENTS[:, 2:end],
        kwargs...
    )
    soilwaterplot, soiltempplot
end

function plot_growth!(plt, model, u, envstart)
    model.dead[] = false
    model.environment_start[] = envstart
    tstop = LIFESPAN * 1.0hr
    prob = DiscreteProblem(model, u, (0.0hr, tstop))
    sol = solve(prob, FunctionMap(scale_by_time = true))
    shootvals = map(i -> sol[ustrip(round(typeof(1hr), i))][:VS] * 25/mol, 1hr:1hr:tstop+1hr)
    rootvals = map(i -> sol[ustrip(round(typeof(1hr), i))][:VR] * -25/mol, 1hr:1hr:tstop+1hr)
    rng = envstart:1.0hr:envstart+tstop
    plot!(plt, rng, shootvals;
          linecolor=:black,
          label="Shoot"
    )
    plot!(plt, rng, rootvals;
        linecolor=:grey,
        label="Root"
    )
    return model, sol
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

" Plot all starting states in all environments"
function plot_years(model, environments, u, envstart)
    locplots = []
    xticks = (collect(0:8760*2:8760*YEARS), string.(collect(STARTYR:2:STOPYR)))
    stateplot = plot(;
        grid=false,
        ylab="Plant structural mass",
        showaxis=false,
        xlims=(-2, -1),
    )
    plot_growth!(stateplot, model, u, envstart)
    microplots = plot_microclim(model, :t1;
         linewidth=MULTIYEARLINEWIDTH,
         xlims=(-100, -10),
         ylabs=(swp="", vpd="", st=""),
         grid=false,
         showaxis=false
    )
    heights, num_plots = calc_heights(microplots)
    legends = plot(stateplot, microplots...;
        legend=:right,
        xticks=xticks,
        link=:x,
        xlab="",
        ylab="",
        margin=0px,
        bottom_margin=20px,
        left_margin=20px,
        layout=grid(num_plots, 1, heights=heights)
    )
    local heights, num_plots
    for (i, envname) in enumerate(keys(environments))
        println(envname)
        model.environment = environments[envname]
        # Allometry Y intercept has to match seed size.
        println("using starting mass given in state for B0")
        if i == 1
            stateplot = plot_months(model, u, envstart; ylab="Structural mass (g)")
        else
            stateplot = plot_months(model, u, envstart; ylab="", yshowaxis=false)
        end
        # xaxis = ((0,8760*11), 0:8760*2:8760*11)
        if i == 1
            microplots = plot_microclim(model, envname; linewidth=MULTIYEARLINEWIDTH, legend=false)
            heights, num_plots = calc_heights(microplots)
            push!(locplots,
                  plot(stateplot, microplots...;
                      xticks=xticks,
                      link=:x,
                      margin=0px,
                      bottom_margin=20px,
                      left_margin=20px,
                      layout=grid(num_plots, 1, heights=heights)
                  )
            )
        else
            microplots = plot_microclim(model, envname; 
                linewidth=MULTIYEARLINEWIDTH,
                ylabs=(swp="", vpd="", st=""), 
                yshowaxis=false, 
                legend=false
            )
            heights, num_plots = calc_heights(microplots)
            push!(locplots,
                plot(stateplot, microplots...;
                    xticks=xticks,
                    link=:x,
                    margin=0px,
                    bottom_margin=20px,
                    left_margin=-80px,
                    layout=grid(num_plots, 1, heights=heights)
                )
            )
        end
    end
    plt = plot(locplots...,
        legends,
        size=(1300,700),
        dpi=DPI,
        layout=grid(1, length(locplots)+1, widths=envwidths(environments))
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
    prob = DiscreteProblem(model, ustrip(u), tspan)
    sol = solve(prob, FunctionMap(scale_by_time = true))
    n = length(u) ÷ 2
    solt = sol'
    solt[:, 4:6] = solt[:, 4:6] * -1
    plts = plot_list(model, solt, envstart, tstop, months, (1,ustrip(tstop)); legend=false)
    legends = plot_list(model, solt, envstart, tstop, months, (-100, -10);
        legend=:bottomright,
        xlims=(-100, -10),
        grid=false,
        showaxis=false,
        ylab=""
    )
    plt = plot(plts, legends;
        layout=grid(1, 2; widths=[0.85, 0.15]),
        size=(1300,700),
        dpi=DPI
    )
    plt
end

function plot_list(model, solt, envstart, tstop, months, xlims; kwargs...)
    rnge = ustrip(envstart:1hr:envstart+tstop)
    yticks = ([10^-1, 10^0, 10^1],
              string.(["-0.1", "-1", "-10"], Ref(" kPa")))
    xticks=(ustrip(0hr:MONTH_HOURS:tstop), MONTHS[1+STARTMONTH:months+1+STARTMONTH])
    rad = radiation(model.environment)
    radplot = plot(rad[rnge];
        ylab="Radiation",
        color=:black,
        kwargs...,
        legend=:none
    )
    solplot = plot(solt;
        linewidth=LINEWIDTH, color=[1 2 3],
        labels=reshape([STATELABELS...], 1, 6),
        ylab="Plant State\nVariables\n(C/N mol)",
        xlab="",
        xshowaxis=false,
        kwargs...
    )
    soilwaterplot, soiltempplot =
        plot_microclim(model, ""; linewidth=MULTIYEARLINEWIDTH, kwargs...)
    tempcorrplot = plot(
        ustrip.(hcat(model.records[1].vars.tempcorrection, model.records[2].vars.tempcorrection));
        ylabel="Temperatuire\ncorrection",
        xlabel="",
        linewidth=LINEWIDTH,
        labels=["Shoot" "Root"],
        xshowaxis=false,
        kwargs...
    )
    rateplot = plot(ustrip.(hcat(model.records[1].vars.rate, model.records[2].vars.rate));
        ylabel="Growth rate\n(mol mol^-1 d^-1)",
        labels=["Shoot" "Root"],
        linewidth=LINEWIDTH,
        xshowaxis=false,
        kwargs...
    )
    depthplot = plot(ustrip.(hcat(model.records[1].vars.height, model.records[2].vars.height .* -1));
        ylabel="Height/\nDepth (m)",
        labels=["Shoot" "Root"],
        linewidth=LINEWIDTH,
        xshowaxis=false,
        kwargs...
    )
    maxswpplot = plot(model.records[1].vars.swp; ylab="Available\nsoil water\npotential",
        xticks=xticks,
        linewidth=LINEWIDTH,
        labels=ENVRANGE,
        xshowaxis=false,
        kwargs...,
        legend=false
    )
    plts = (solplot, depthplot, rateplot, soiltempplot, tempcorrplot, maxswpplot, soilwaterplot, radplot)
    x = 0.8/8
    plot(plts...;
        xlims=xlims,
        xticks=xticks,
        link=:x,
        layout=grid(length(plts), 1, heights=[0.2, x, x, x+0.03, x-0.03, x+0.03, x-0.03, x, x])
    )
end

function plot_assim(model, environment, u, envstart, months, plotids)
    tstop = round(typeof(1hr), months * MONTH_HOURS)
    rnge = ustrip(envstart:1hr:envstart+tstop)
    plotstart = 1hr
    model.environment = environment
    model.environment_start[] = envstart
    model.dead[] = false
    tspan = ustrip.((plotstart, tstop))
    prob = DiscreteProblem(model, ustrip.(u), tspan)
    sol = solve(prob, FunctionMap(scale_by_time = true))
    xticks=(ustrip(0hr:MONTH_HOURS:tstop), MONTHS[1+STARTMONTH:months+1+STARTMONTH])
    yticks = (-0.05:0.05:0.2)
    st = soiltemperature(model.environment) .|> °C
    rh = relhumidity(model.environment)
    airtemp = airtemperature(model.environment)
    vpd = uconvert.(u"kPa", vapour_pressure_deficit.(airtemp, rh))
    n = length(u) ÷ 2
    solt = sol' .* 25g
    solt[:, 4:6] = solt[:, 4:6] * -1
    solplot = plot(solt,
        legend=:topleft,
        linewidth=LINEWIDTH,
        color=[1 2 3],
        labels=reshape([STATELABELS...], 1, 6),
        ylabel="State",
        ytick=yticks,
        xticks=xticks,
        xlabel=""
    )
    assimplot = plot(ustrip.(hcat(model.records[1].J[2,1,:]));
        ylabel="Assimilation\n(C-mol/hr)",
        xticks=xticks,
        linewidth=LINEWIDTH,
        legend=false
    )
    avswpplot = plot(hcat(model.records[1].vars.swp);
        ylabel="Available\nsoil water\npotential",
        xticks=xticks,
        linewidth=LINEWIDTH,
        labels=ENVRANGE,
        legend=false
    )
    # radplot = plot(radiation(environment); ylab="Radiation", color=:black, legend=:none, xlabel="Month")
    soilcolors = permutedims(reverse(get.(Ref(ColorSchemes.copper), 0.0:1/7:1)))
    swpyticks = ([10^0, 10^1], string.(["-10e0", "-10e1"], Ref(" kPa")))
    swp = -1 .* ustrip(soilwaterpotential(model.environment)) # .|> MPa
    depthplot = plot(ustrip.(hcat(model.records[1].vars.height, model.records[2].vars.height .* -1));
        ylabel="Height/\nDepth (m)",
        labels=["Shoot" "Root"],
        linewidth=LINEWIDTH,
        xshowaxis=false,
        legend=:topleft
    )
    soilwaterplot = plot(swp[rnge, 2:end];
        yscale=:log10,
        ylab=ENVLABELS[:swp],
        xticks=xticks,
        yflip=true,
        yticks=swpyticks,
        color=soilcolors,
        labels=ENVINCREMENTS[:, 2:end]
    )
    plts = (solplot, assimplot, avswpplot, depthplot, soilwaterplot)[plotids]
    heights = assimheights(plts, plotids)
    plt = plot(plts...;
        xlims=(1,ustrip(tstop)),
        left_margin=60px,
        link=:x,
        layout=grid(length(plts), 1, heights=heights),
        size=(1300,700),
        dpi=DPI
    )
    plt
end

function assimheights(plts, rnge)
    l = length(plts)
    if rnge.start == 1
        h = 1/(l+1)
        [2h, [h for x in 1:l-1]...]
    else
        [1/l for x in plts]
    end
end

function plot_scaling(model, environment, u, envstart, months, x)
    tstop = round(typeof(1hr), months * MONTH_HOURS)
    plotstart = 1hr
    model.environment = environment
    model.environment_start[] = envstart
    model.dead[] = false
    tstop = round(typeof(1hr), months * MONTH_HOURS)
    tspan = ustrip.((plotstart, tstop))
    scaling = oneunit(model.params[2].scaling_pars.M_Vscaling)
    xticks=(ustrip(0hr:MONTH_HOURS:tstop), MONTHS[1+STARTMONTH:months+1+STARTMONTH])
    plt = plot(
        legend=:topleft,
        linewidth=LINEWIDTH,
        labels=reshape([STATELABELS...], 1, 6),
        ylabel="State\n(C/N mol)",
        xlabel=""
    )
    rnge = collect(scaling/x:(scaling-scaling/x)/4:scaling)
    for (i, s) in enumerate(rnge)
        m = @set model.params[2].scaling_pars.M_Vscaling = s
        prob = DiscreteProblem(m, ustrip.(u), tspan)
        sol = solve(prob, FunctionMap(scale_by_time = true))
        n = length(u) ÷ 2
        solt = sol'
        solt[:, 4:6] = solt[:, 4:6] * -1
        alpha = (i/length(rnge))
        println(alpha)
        plot!(plt, solt, color=[1 2 3], alpha=ALPHA)
        annotate!(plt, ustrip(tstop), solt[end, 1] - 0.05, text(string(round(ustrip(s), digits=2)), 7))
    end
    plot(plt; xlims=(1, ustrip(tstop)), size=(1300,700), dpi=DPI, legend=:none)
end
