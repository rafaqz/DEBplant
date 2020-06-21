dir = "DEBSCRIPTS" in keys(ENV) ? ENV["DEBSCRIPTS"] : pwd()
include(joinpath(dir, "load.jl"))
include(joinpath(dir, "plantstates.jl"))


using PlotThemes, LaTeXStrings, Dates

MONTHS = [Dates.LOCALES["english"].months... Dates.LOCALES["english"].months...]
MONTH_HOURS = 365.25 / 12 * 24hr
STARTMONTH = 8
YLIMS = (-3,4)
LIFESPAN = 4380
LINEWIDTH = 1.3
MULTIYEARLINEWIDTH = 0.2
ENVLABS = ("Soil water\npotential", "Vapour press.\ndeficit", "Soil\ntemperature")
DPI = 150 
STARTYR = 2005 
STOPYR = 2011 
YEARS = STOPYR-STARTYR
SMALLFONT = 6
ALPHA = 1.0


import Base: round
round(::Type{T}, x::Quantity) where {T<:Quantity} = T(round(typeof(one(T)), uconvert(unit(T), x).val))

function plot_microclim(model, envname, rnge = 1:1:size(radiation(model.environment), 1);
                        ylabs=ENVLABS, kwargs...)
    swp = -1 .* ustrip(soilwaterpotential(model.environment)) # .|> MPa
    st = soiltemperature(model.environment) .|> °C
    vpd = uconvert.(u"kPa", vapour_pressure_deficit.(airtemperature(model.environment), relhumidity(model.environment)))
    soilcolors = permutedims(reverse(get.(Ref(ColorSchemes.copper), 0.0:1/7:1)))
    yticks = ([10^-1, 10^0, 10^1, 10^2, 10^3, 10^4],
              string.(["-10e-1", "-10e0", "-10e1", "-10e2", "-10e3", "-10e4"], Ref(" kPa")))
    soilwaterplot = plot(rnge * hr, swp[rnge, 2:end]; yscale=:log10,
                         ylab=ylabs[1], yflip=true, yticks=yticks,
                         color=soilcolors, labels=INCREMENTS[:, 2:end], kwargs...)
    # Make the plot full range keeping log scale. But not for the legend...
    if !(:showaxis in keys(kwargs))
        plot!(soilwaterplot, (1, 20000.0))
        plot!(soilwaterplot, (1, 0.1))
    end
    soiltempplot = plot(rnge * hr, st[rnge, 2:end];  ylab=ylabs[3], ylims=(-2, 50),
                        color=soilcolors, xlab=envname, labels=INCREMENTS[:, 2:end], kwargs...)
    soilwaterplot, soiltempplot
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
function plot_years(model, environments, states, envstart, plotname)
    locplots = []
    statename = first(keys(states))
    xticks = (collect(0:8760*2:8760*YEARS), string.(collect(STARTYR:2:STOPYR)))

    stateplot = plot(; grid=false, showaxis=false, xlims=(-2, -1), legend=:bottomright)
    plot_sol!(stateplot, model, states[statename], envstart)
    microplots = plot_microclim(model, :t1; linewidth=MULTIYEARLINEWIDTH,
                                xlims=(-100, -10), ylabs=("", "", ""), grid=false,
                                showaxis=false)
    heights, num_plots = calc_heights(microplots)
    legends = plot(stateplot, microplots...;
                  xticks=xticks, link=:x, margin=0px, bottom_margin=20px, left_margin=20px,
                  layout=grid(num_plots, 1, heights=heights))

    local heights, num_plots
    for (i, envname) in enumerate(keys(environments))
        println(envname)
        model.environment = environments[envname]
        # Allometry Y intercept has to match seed size.
        println("using starting mass given in state for B0")
        if i == 1
            stateplot = plot_months("Shoot and Root\nstructural mass", model, states[statename], envstart)
        else
            stateplot = plot_months("", model, states[statename], envstart; yshowaxis=false)
        end
        # xaxis = ((0,8760*11), 0:8760*2:8760*11)
        if i == 1
            microplots = plot_microclim(model, envname; linewidth=MULTIYEARLINEWIDTH, legend=false)
            heights, num_plots = calc_heights(microplots)
            push!(locplots, plot(stateplot, microplots...;
                  xticks=xticks, link=:x, margin=0px, bottom_margin=20px, left_margin=20px,
                  layout=grid(num_plots, 1, heights=heights)))
        else
            microplots = plot_microclim(model, envname; linewidth=MULTIYEARLINEWIDTH,
                                        ylabs=("", "", ""), yshowaxis=false, legend=false)
            heights, num_plots = calc_heights(microplots)
            push!(locplots, plot(stateplot, microplots...;
                 xticks=xticks, link=:x, margin=0px, bottom_margin=20px, left_margin=-80px,
                 layout=grid(num_plots, 1, heights=heights)))
        end

    end

    plt = plot(locplots..., legends, size=(1300,700), dpi=DPI,
               layout=grid(1, length(locplots)+1, widths=envwidths(environments)))
    savefig("plots/$plotname")
    plt
end

envwidths(envs) = begin
    n = length(keys(envs))
    [[0.9/n for i in 1:n]..., 0.1]
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
    solt = sol'
    solt[:, 4:6] = solt[:, 4:6] * -1
    plts = plot_list(model, solt, envstart, tstop, months, (1,ustrip(tstop)); legend=false)
    legends = plot_list(model, solt, envstart, tstop, months, (-100, -10); legend=:bottomright,
                        xlims=(-100, -10), grid=false, showaxis=false, ylab="")
    plt = plot(plts, legends, layout=grid(1, 2; widths=[0.85, 0.15]), size=(1300,700), dpi=DPI)
    savefig("plots/crossover")
    plt
end

plot_list(model, solt, envstart, tstop, months, xlims; kwargs...) = begin
    rnge = ustrip(envstart:1hr:envstart+tstop)
    yticks = ([10^-1, 10^0, 10^1],
              string.(["-0.1", "-1", "-10"], Ref(" kPa")))
    xticks=(ustrip(0hr:MONTH_HOURS:tstop), MONTHS[1+STARTMONTH:months+1+STARTMONTH])
    rad = radiation(model.environment)
    radplot = plot(rad[rnge]; ylab="Radiation", color=:black, kwargs..., legend=:none)
    solplot = plot(solt; linewidth=LINEWIDTH, color=[1 2 3],
                   labels=reshape([STATELABELS...], 1, 6), ylab="Plant State\nVariables\n(C/N mol)",
                   xlab="", xshowaxis=false, kwargs...)
    soilwaterplot, soiltempplot = 
        plot_microclim(model, ""; linewidth=MULTIYEARLINEWIDTH, kwargs...)
    tempcorrplot = plot(ustrip.(hcat(model.records[1].vars.tempcorrection, model.records[2].vars.tempcorrection));
                    ylabel="Temperatuire\ncorrection", xlabel="", linewidth=LINEWIDTH,
                    labels=["Shoot" "Root"], xshowaxis=false, kwargs...)
    rateplot = plot(ustrip.(hcat(model.records[1].vars.rate, model.records[2].vars.rate));
                    ylabel="Growth rate\n(mol mol^-1 d^-1)",
                    labels=["Shoot" "Root"], linewidth=LINEWIDTH, xshowaxis=false, kwargs...)
    depthplot = plot(ustrip.(hcat(model.records[1].vars.height, model.records[2].vars.height .* -1));
                     ylabel="Height/\nDepth (m)", labels=["Shoot" "Root"], 
                     linewidth=LINEWIDTH, xshowaxis=false, kwargs...)
    maxswpplot = plot(model.records[1].vars.swp; ylab="Available\nsoil water\npotential",
                      xticks=xticks, linewidth=LINEWIDTH, labels=RANGE, 
                      xshowaxis=false, kwargs..., legend=false)

    plts = (solplot, depthplot, rateplot, soiltempplot, tempcorrplot, soilwaterplot, maxswpplot, radplot)
    x = 0.8/8
    plot(plts...; xlims=xlims, xticks=xticks, link=:x,
         layout=grid(length(plts), 1, heights=[0.2, x, x, x+0.03, x-0.03, x+0.03, x-0.03, x, x]))
end

function plot_assim(model, environment, u, envstart, months, plotids, name)
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
    # yticks = (-0.05:0.05:0.2)
    st = soiltemperature(model.environment) .|> °C
    vpd = uconvert.(u"kPa", vapour_pressure_deficit.(airtemperature(model.environment), relhumidity(model.environment)))
    n = length(u) ÷ 2
    solt = sol' .* 25g
    solt[:, 4:6] = solt[:, 4:6] * -1
    solplot = plot(solt, legend=:topleft, linewidth=LINEWIDTH, color=[1 2 3],
                   labels=reshape([STATELABELS...], 1, 6), ylabel="State", ytick=yticks,
                   xticks=xticks, xlabel="")
    assimplot = plot(ustrip.(hcat(model.records[1].J[2,1,:])); 
                     ylabel="Assimilation\n(C-mol/hr)", xticks=xticks, linewidth=LINEWIDTH, legend=false)
    avswpplot = plot(hcat(model.records[1].vars.swp); ylabel="Available\nsoil water\npotential",
                   xticks=xticks, linewidth=LINEWIDTH,
                   labels=RANGE, legend=false)
    # radplot = plot(radiation(environment); ylab="Radiation", color=:black, legend=:none, xlabel="Month")
    soilcolors = permutedims(reverse(get.(Ref(ColorSchemes.copper), 0.0:1/7:1)))
    swpyticks = ([10^0, 10^1],
              string.(["-10e0", "-10e1"], Ref(" kPa")))
    swp = -1 .* ustrip(soilwaterpotential(model.environment)) # .|> MPa
    depthplot = plot(ustrip.(hcat(model.records[1].vars.height, model.records[2].vars.height .* -1));
                     ylabel="Height/\nDepth (m)", labels=["Shoot" "Root"], 
                     linewidth=LINEWIDTH, xshowaxis=false, legend=:topleft)
    soilwaterplot = plot(swp[rnge, 2:end]; yscale=:log10, 
                         ylab=ENVLABS[1], xticks=xticks, yflip=true, yticks=swpyticks,
                         color=soilcolors, labels=INCREMENTS[:, 2:end]) 
    plts = (solplot, assimplot, avswpplot, depthplot, soilwaterplot)[plotids]
    heights = assimheights(plts, plotids)
    plt = plot(plts...; xlims=(1,ustrip(tstop)), left_margin=60px, link=:x,
               layout=grid(length(plts), 1, heights=heights), 
               size=(1300,700), dpi=DPI)
    savefig("plots/$name")
    plt
end

assimheights(plts, rnge) = begin
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
    plt = plot(legend=:topleft, linewidth=LINEWIDTH, labels=reshape([STATELABELS...], 1, 6),
               ylabel="State\n(C/N mol)", xlabel="")
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
model = set_allometry(models[:bb], u)
INCREMENTS = reshape([string.(Microclimate.get_increments(environment))...], 1, 8)
RANGE = reshape([string.(Microclimate.get_range(environment))...], 1, 2)

temps = [(0:50)...] * °C
plot(x -> tempcorr(tempcorr_pars(model.shared), K(x)), temps,
     legend=false, ylabel="Correction", xlabel="Temperature")
savefig("plots/tempcorr")


plot_crossover(model, environment, u, envstart, 2)
plot_assim(model, environment, u, envstart, 6, 1:1, "growth")
plot_assim(model, environment, u, envstart, 6, 1:2, "growthassim")
plot_assim(model, environment, u, envstart, 6, 1:3, "assimall")
plot_assim(model, environment, u, envstart, 1, 1:1, "assimshort")
plot_assim(model, environment, u, envstart, 6, 3:5, "assimswp")
plot_assim(model, environment, u, envstart, 6, 1:5, "assimall")

model = set_allometry(models[:bb], u)
scaling = plot_scaling(model, environment, u, envstart, 12, 1.5); 
gr() 
envstart = 1hr 
crossover
plot_years(model, environments, states, envstart, "all") 
plot_years(model, Dict(:t1=>environments[:t1]), states, envstart, "t1") 
plot_years(model, Dict(:t2=>environments[:t2]), states, envstart, "t2") 
plot_years(model, Dict(:t3=>environments[:t3]), states, envstart, "t3") 

using Microclimate
@set! model.environment = MicroclimControl(1000.0W*m^-2, 0.0m, 290.0K, 0.8, 0.5m*s^-1, 290.0K, -100.0kPa, 0.7)
assim = plot_assim(model, MicroclimControl(), u, envstart, 6, 1:1, "constantgrowth")
