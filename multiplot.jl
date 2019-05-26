dir = "DEBSCRIPTS" in keys(ENV) ? ENV["DEBSCRIPTS"] : pwd()
include(joinpath(dir, "load.jl"))
include(joinpath(dir, "plantstates.jl"))


using PlotThemes, LaTeXStrings, Dates

MONTHS = [Dates.LOCALES["english"].months... Dates.LOCALES["english"].months...]
MONTH_HOURS = 365.25 / 12 * 24hr
STARTMONTH = 8
YLIMS = (-3,7)
LIFESPAN = 4380
LINEWIDTH = 1.3
MULTIYEARLINEWIDTH = 0.2
ENVLABS = ("Soil water\npotential", "Vapour press.\ndeficit", "Soil\ntemperature")
DPI = 150
STARTYR = 2005
STOPYR = 2011
YEARS = STOPYR-STARTYR
SMALLFONT = 6
alpha = 0.9


import Base: round
round(::Type{T}, x::Quantity) where {T<:Quantity} = T(round(typeof(one(T)), uconvert(unit(T), x).val))

function plot_microclim(model, envname, rnge = 1:1:size(radiation(model.environment), 1);
                        ylabs=ENVLABS, kwargs...)
    swp = -1 .* ustrip(soilwaterpotential(model.environment)) # .|> MPa
    st = soiltemperature(model.environment) .|> °C
    vpd = uconvert.(u"kPa", vapour_pressure_deficit.(airtemperature(model.environment), relhumidity(model.environment)))
    soilcolors = permutedims(reverse(get.(Ref(ColorSchemes.copper), 0.0:1/7:1)))
    # yticks = ["-10⁰ KPa", "-10¹ KPa", "-10² KPa", "-10³ KPa", "-10⁴ KPa"]
    # yticks = [-10^0u"kPa", -10^1u"kPa", -10^2u"kPa", -10^3u"kPa", -10^4u"kPa"]
    yticks = ([10^-1, 10^0, 10^1, 10^2, 10^3, 10^4],
              string.(["-0.1", "-1", "-10", "-10^2", "-10^3", "-10^4"], Ref(" kPa")))
    soilwaterplot = plot(rnge * hr, swp[rnge, 2:end]; yscale=:log10,
                         ylab=ylabs[1], yflip=true, yticks=yticks,
                         color=soilcolors, labels=INCREMENTS[:, 2:end], kwargs...)
    # Make the plot full range keeping log scale. But not for the legend...
    if !(:showaxis in keys(kwargs))
        plot!(soilwaterplot, (1, 20000.0))
        plot!(soilwaterplot, (1, 0.1))
    end
    # plot!(soilwaterplot, (-100, 0))
    soiltempplot = plot(rnge * hr, st[rnge, 2:end];  ylab=ylabs[3], ylims=(-2, 50),
                        color=soilcolors, labels=INCREMENTS[:, 2:end], kwargs...)
    vpdplot = plot(rnge * hr, vpd[rnge, :]; ylab=ylabs[2], alpha=alpha,
                   xlab=titlecase(string(envname)), ylims=(-550, 0),
                  labels=RANGE, color=[:black :grey], kwargs...)
    soilwaterplot, soiltempplot, vpdplot
end

function plot_sol!(plt, model, u, envstart)
    model.dead[] = false
    model.environment_start[] = envstart
    tstop = LIFESPAN * hr
    prob = DiscreteProblem(model, u, (0hr, tstop))
    sol = solve(prob, FunctionMap(scale_by_time = true))
    shootvals = map(i -> sol[ustrip(round(typeof(1hr), i))][:VS] * 25g/mol, 1hr:1hr:tstop+1hr)
    rootvals = map(i -> sol[ustrip(round(typeof(1hr), i))][:VR] * -25g/mol, 1hr:1hr:tstop+1hr)
    rng = envstart:1hr:envstart+tstop
    plot!(plt, rng, shootvals; linecolor=:black, label="Shoot")
    plot!(plt, rng, rootvals; linecolor=:grey, label="Root")
    return model, sol
end

" Plot growth starting every moth in the available timespan"
function plot_months(title, model, u, envstart; kwargs...)
    solplots = plot(; ylab=title, ylims=YLIMS,
                    legend=:none, link=:x, kwargs...)
    for i in 1:YEARS*12
        println("month: ", i)
        envstart += MONTH_HOURS
        # Round to the start of a day
        envstart_hour = round(typeof(1hr), round(typeof(1d), envstart)) + 1hr
        plot_sol!(solplots, model, u, envstart_hour)
    end
    solplots
end

" Plot all starting states in all environments"
function plot_years(model, environments, states, envstart)
    locplots = []
    statename = first(keys(states))
    xticks = (collect(0:8760*2:8760*YEARS), string.(collect(STARTYR:2:STOPYR)))

    stateplot = plot(; grid=false, showaxis=false, xlims=(-2, -1))
    plot_sol!(stateplot, model, states[statename], envstart)
    microplots = plot_microclim(model, :t1; linewidth=MULTIYEARLINEWIDTH,
                                xlims=(-100, -10), ylabs=("", "", ""), grid=false,
                                showaxis=false)
    heights, num_plots = calc_heights(microplots)
    legends = plot(stateplot, microplots...;
                  xticks=xticks, link=:x, margin=0px, left_margin=20px,
                  layout=grid(num_plots, 1, heights=heights))

    local heights, num_plots
    for (i, envname) in enumerate(keys(environments))
        println(envname)
        model.environment = environments[envname]
        # Allometry Y intercept has to match seed size.
        println("using starting mass given in state for B0")
        if i == 1
            stateplot = plot_months("Shoot and Root\nmass", model, states[statename], envstart)
        else
            stateplot = plot_months("", model, states[statename], envstart; yshowaxis=false)
        end
        # xaxis = ((0,8760*11), 0:8760*2:8760*11)
        if i == 1
            microplots = plot_microclim(model, envname; linewidth=MULTIYEARLINEWIDTH, legend=false)
            heights, num_plots = calc_heights(microplots)
            push!(locplots, plot(stateplot, microplots...;
                  xticks=xticks, link=:x, margin=0px, right_margin=0px,
                  layout=grid(num_plots, 1, heights=heights)))
        else
            microplots = plot_microclim(model, envname; linewidth=MULTIYEARLINEWIDTH,
                                        ylabs=("", "", ""), yshowaxis=false, legend=false)
            heights, num_plots = calc_heights(microplots)
            push!(locplots, plot(stateplot, microplots...;
                 xticks=xticks, link=:x, margin=0px, left_margin=-80px,
                 layout=grid(num_plots, 1, heights=heights)))
        end

    end

    plt = plot(locplots..., legends, size=(1000,1000), dpi=DPI,
               layout=grid(1, length(locplots)+1, widths=[0.3, 0.3, 0.3, 0.1]))
    savefig("plots/all")
    plt
end

calc_heights(microplots) = begin
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
    rnge = ustrip(envstart:1hr:envstart+tstop)
    xticks=(ustrip(0hr:MONTH_HOURS:tstop), MONTHS[1+STARTMONTH:months+1+STARTMONTH])
    swp = -1 .* ustrip(soilwaterpotential(model.environment)) # .|> MPa
    st = soiltemperature(model.environment) .|> °C
    vpd = vapour_pressure_deficit.(airtemperature(model.environment), relhumidity(model.environment))
    rad = radiation(model.environment)
    soilcolors = permutedims(reverse(get.(Ref(ColorSchemes.copper), 0.0:1/7:1)))
    solt = sol'
    solt[:, 4:6] = solt[:, 4:6] * -1
    solplot = plot(solt, linewidth=LINEWIDTH, color=[1 2 3],
                   labels=reshape([STATELABELS...], 1, 6), ylabel="State\n(C/N mol)",
                   xlabel="", xshowaxis=false)
    radplot = plot(rad[rnge]; ylab="Radiation", color=:black, legend=:none, xshowaxis=false)
    yticks = ([10^-1, 10^0, 10^1],
              string.(["-0.1", "-1", "-10"], Ref(" kPa")))
    soilwaterplot = plot(swp[rnge, 2:end]; yflip=true, yticks=yticks, yscale=:log10, ylab=ENVLABS[1],
                         color=soilcolors, labels=INCREMENTS[:, 2:end], legendfontsize=SMALLFONT,
                         linewidth=LINEWIDTH, alpha=alpha, xshowaxis=false)
    soiltempplot = plot(st[rnge, 2:end]; ylab=ENVLABS[3], legendfontsize=SMALLFONT,
                        color=soilcolors, linewidth=LINEWIDTH, xshowaxis=false,
                        labels=INCREMENTS[:, 2:end], alpha=alpha)
    vpdplot = plot(uconvert.(u"kPa", vpd[rnge, :]); ylab=ENVLABS[2], alpha=alpha,
                  labels=RANGE, linewidth=LINEWIDTH, xlabel="Month", legendfontsize=SMALLFONT)
    tempcorrplot = plot(ustrip.(hcat(model.records[1].vars.tempcorrection, model.records[2].vars.tempcorrection));
                    ylabel="Temperatuire\ncorrection", xlabel="", linewidth=LINEWIDTH,
                    alpha=alpha, legendfontsize=SMALLFONT,
                    labels=["Shoot" "Root"], xshowaxis=false)
    rateplot = plot(ustrip.(hcat(model.records[1].vars.rate, model.records[2].vars.rate));
                    ylabel="Growth rate\n(mol mol^-1 d^-1)",
                    labels=["Shoot" "Root"], linewidth=LINEWIDTH, legendfontsize=SMALLFONT,
                    alpha=alpha, xshowaxis=false)
    depthplot = plot(ustrip.(hcat(model.records[1].vars.height, model.records[2].vars.height .* -1));
                     ylabel="Height/\nDepth (m)", labels=["Shoot" "Root"], legendfontsize=SMALLFONT,
                     linewidth=LINEWIDTH, alpha=alpha, xshowaxis=false)
    maxswpplot = plot(model.records[1].vars.swp; yaxis=("Available\nsoil water\npotential"),
                      xticks=xticks, linewidth=LINEWIDTH, labels=RANGE, legend=false, xshowaxis=false)

    plts = (solplot, depthplot, rateplot, soiltempplot, tempcorrplot, soilwaterplot, maxswpplot, radplot, vpdplot)
    x = 0.8/8
    plt = plot(plts...; xlims=(1,ustrip(tstop)), xticks=xticks, link=:x,
               layout=grid(length(plts), 1, heights=[0.2, x, x, x+0.03, x-0.03, x+0.03, x-0.03, x, x]),
               size=(1000,1200), dpi=DPI)
    savefig("plots/crossover")
    plt
end

function plot_assim(model, environment, u, envstart, months)
    tstop = round(typeof(1hr), months * MONTH_HOURS)
    plotstart = 1hr
    model.environment = environment
    model.environment_start[] = envstart
    model.dead[] = false
    tspan = ustrip.((plotstart, tstop))
    prob = DiscreteProblem(model, ustrip(u), tspan)
    sol = solve(prob, FunctionMap(scale_by_time = true))
    xticks=(ustrip(0hr:MONTH_HOURS:tstop), MONTHS[1+STARTMONTH:months+1+STARTMONTH])
    n = length(u) ÷ 2
    solt = sol'
    solt[:, 4:6] = solt[:, 4:6] * -1
    solplot = plot(solt, legend=:topleft, linewidth=LINEWIDTH, color=[1 2 3],
                   labels=reshape([STATELABELS...], 1, 6), ylabel="State\n(C/N mol)",
                   xticks=xticks, xlabel="")
    assimplot = plot(ustrip.(hcat(model.records[1].J[2,1,:]));
                     ylabel="Assimilation\n(C-mol hr^-1)", xticks=xticks, linewidth=LINEWIDTH, legend=false)
    swpplot = plot(hcat(model.records[1].vars.swp); yaxis=("Available\nsoil water\npotential"),
                   xlabel="Month", xticks=xticks, linewidth=LINEWIDTH,
                   labels=RANGE, legend=false)
    radplot = plot(radiation(environment); ylab="Radiation", color=:black, legend=:none, xlabel="Month")
    plt = plot(solplot, assimplot, swpplot; xlims=(1,ustrip(tstop)), link=:x,
         layout=grid(3, 1, heights=[0.5, 0.25, 0.25]), size=(1000,1000), dpi=DPI)
    savefig("plots/assim")
    plt
end

function plot_scaling(model, environment, u, envstart, months, x)
    tstop = round(typeof(1hr), months * MONTH_HOURS)
    plotstart = 1hr
    model.environment = environment
    model.environment_start[] = envstart
    model.dead[] = false
    tstop = round(typeof(1hr), months * MONTH_HOURS)
    tspan = ustrip.((plotstart, tstop))
    scaling = oneunit(model.params[2].shape_pars.M_Vscaling)
    xticks=(ustrip(0hr:MONTH_HOURS:tstop), MONTHS[1+STARTMONTH:months+1+STARTMONTH])
    plt = plot(legend=:topleft, linewidth=LINEWIDTH, labels=reshape([STATELABELS...], 1, 6),
               ylabel="State\n(C/N mol)", xlabel="")
    rnge = collect(scaling/x:(scaling-scaling/x)/4:scaling)
    for (i, s) in enumerate(rnge)
        m = @set model.params[2].shape_pars.M_Vscaling = s
        prob = DiscreteProblem(m, ustrip.(u), tspan)
        sol = solve(prob, FunctionMap(scale_by_time = true))
        n = length(u) ÷ 2
        solt = sol'
        solt[:, 4:6] = solt[:, 4:6] * -1
        alpha = (i/length(rnge))
        println(alpha)
        plot!(plt, solt, color=[1 2 3], alpha=alpha)
        annotate!(plt, ustrip(tstop), solt[end, 1] - 0.05, text(string(round(ustrip(s), digits=2)), 7))
    end
    plot(plt; xlims=(1, ustrip(tstop)), size=(900,900), dpi=DPI, legend=:none)
    savefig("plots/scaling")
    plt
end

theme(:sand)
theme(:solarized)
theme(:juno)
theme(:solarized_light)
theme(:wong2)
x = 1.5

# pyplot()
# pyplot()
# plotly()

environments, _ = loadenvironments(dir)
models = OrderedDict()
modeldir = joinpath(dir, "models")
include.(joinpath.(Ref(modeldir), readdir(modeldir)));
envstart = round(typeof(1hr), STARTMONTH * MONTH_HOURS)
u = states["Large seed"]
# model = models[:bb]
environment = environments[:t1]
model = set_allometry(models[:bbiso], u)
INCREMENTS = reshape([string.(Microclimate.get_increments(environment))...], 1, 8)
RANGE = reshape([string.(Microclimate.get_range(environment))...], 1, 2)

crossover = plot_crossover(model, environment, u, envstart, 2);
crossover
assim = plot_assim(model, environment, u, envstart, 6);

model = set_allometry(models[:bb], u)
scaling = plot_scaling(model, environment, u, envstart, 12, 1.5);

gr()
envstart = 1hr
allplot = plot_years(model, environments, states, envstart)
allplot
