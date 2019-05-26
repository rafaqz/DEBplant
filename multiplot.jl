dir = "DEBSCRIPTS" in keys(ENV) ? ENV["DEBSCRIPTS"] : pwd()
include(joinpath(dir, "load.jl"))
include(joinpath(dir, "plantstates.jl"))


using PlotThemes, LaTeXStrings, Dates

MONTHS = Dates.LOCALES["english"].months
MONTH_HOURS = 365.25 / 12 * 24hr
YLIMS = (-3,7)
LIFESPAN = 4380
LINEWIDTH = 1.3
MULTIYEARLINEWIDTH = 0.2
ENVLABS = ("Soil water\npotential", "Vapour press.\ndeficit", "Soil\ntemperature")
DPI = 150
STARTYR = 2005
STOPYR = 2011
YEARS = STOPYR-STARTYR
alpha = 0.7


import Base: round
round(::Type{T}, x::Quantity) where {T<:Quantity} = T(round(typeof(one(T)), uconvert(unit(T), x).val)) 

function plot_microclim(model, envname, rnge = 1:1:size(radiation(model.environment), 1); 
                        ylabs=ENVLABS, kwargs...)
    swp = -1 .* ustrip(soilwaterpotential(model.environment)) # .|> MPa
    st = soiltemperature(model.environment) .|> °C
    vpd = vapour_pressure_deficit.(airtemperature(model.environment), relhumidity(model.environment))
    soilcolors = permutedims(reverse(get.(Ref(ColorSchemes.copper), 0.0:1/7:1)))
    # yticks = ["-10⁰ KPa", "-10¹ KPa", "-10² KPa", "-10³ KPa", "-10⁴ KPa"]
    # yticks = [-10^0u"kPa", -10^1u"kPa", -10^2u"kPa", -10^3u"kPa", -10^4u"kPa"]
    yticks = ([10^0, 10^1, 10^2, 10^3, 10^4], 
              string.(Ref("-10"), ["⁰", "¹", "²", "³", "⁴"], Ref(" kPa")))
    soilwaterplot = plot(rnge * hr, swp[rnge, :]; yflip=true, yscale=:log10, ylab=ylabs[1], 
                         color=soilcolors, yticks=yticks, 
                         labels=INCREMENTS, kwargs...)
    soiltempplot = plot(rnge * hr, st[rnge, :];  ylab=ylabs[3], ylims=(-10, 80), color=soilcolors, 
                        labels=INCREMENTS, kwargs..., legend=false)
    vpdplot = plot(rnge * hr, vpd[rnge, :]; ylab=ylabs[2], alpha=alpha, xlab=titlecase(string(envname)),
                  labels=RANGE, color=[:black :darkgrey], kwargs...)
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
    plot!(plt, rng, shootvals, linecolor=:black)
    plot!(plt, rng, rootvals, linecolor=:grey)
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
    for (i, envname) in enumerate(keys(environments))
        println(envname)
        stateplots = []
        for statename in keys(states)
            println("Plotting: ", statename)
            println(statename)
            model.environment = environments[envname]
            # Allometry Y intercept has to match seed size.
            println("using starting mass given in state for B0")
            if i == 1
                push!(stateplots, plot_months("Shoot and Root\nmass", model, states[statename], envstart))
            else
                push!(stateplots, plot_months("", model, states[statename], envstart; yshowaxis=false))
            end
        end
        # xaxis = ((0,8760*11), 0:8760*2:8760*11)
        xticks = (collect(0:8760*2:8760*YEARS), string.(collect(STARTYR:2:STOPYR)))
        if i == 1
            microplots = plot_microclim(model, envname; linewidth=MULTIYEARLINEWIDTH, 
                                        legend=:topleft, legendfontsize=8, background_legend=RGBA(1,1,1,0.5))
            heights, num_plots = calc_heights(stateplots, microplots)
            push!(locplots, plot(stateplots..., microplots...;
                 xticks=xticks, link=:x, margin=0px, right_margin=0px,
                 layout=grid(num_plots, 1, heights=heights)))
        else
            microplots = plot_microclim(model, envname; linewidth=MULTIYEARLINEWIDTH, 
                                        ylabs=("", "", ""), yshowaxis=false, legend=false)
            heights, num_plots = calc_heights(stateplots, microplots)
            push!(locplots, plot(stateplots..., microplots...;
                 xticks=xticks, link=:x, margin=0px, left_margin=-80px,
                 layout=grid(num_plots, 1, heights=heights)))
        end
    end
    println("Plotting output...")
    plt = plot(locplots..., size=(1200,1200), dpi=DPI, layout=grid(1, length(locplots)))
    savefig("plots/all")
    plt
end

calc_heights(stateplots, microplots) = begin
    num_states = length(stateplots)
    num_microplots = length(microplots)
    num_plots = num_states + num_microplots
    h = 1/(num_microplots + 2num_states)
    heights = vcat([2h for x = 1:num_states], [h for x in 1:num_microplots])
    num_plots = num_states + num_microplots
    heights, num_plots
end

function plot_single(model, environments, states, envstart)
    name = "Shoot and Root mass"
    start = round(Int,ustrip(envstart))
    microplots = plot_microclim(model, "T1", start:1:LIFESPAN+start; margin=0px)
    model.environment = environments[:t2]
    solplot = plot(yaxis=(name), legend=:none, link=:x)
    model, sol = plot_sol!(solplot, model, states[name], envstart)
    rateplot = plot(hcat(model.records[1].vars.rate, model.records[2].vars.rate); 
                    yaxis=("Growth rate"), labels=["Shoot" "Root"])
    assimplot = plot(model.records[1].J[4,1,:]; yaxis=("Assimilation"))
    plt = plot(solplot, rateplot, assimplot, microplots...; xaxis=(0:MONTH_HOURS:LIFESPAN, MONTHS), 
               link=:x, layout=grid(6, 1, heights=[0.3, 0.2, 0.15, 0.15, 0.1, 0.1]), size=(600,600), dpi=DPI)
    savefig("plots/single")
    plt
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
    xticks=(ustrip(0hr:round(typeof(1hr), MONTH_HOURS):tstop), MONTHS[1:months+1])
    solplot = plot(solt, legend=:topright, linewidth=LINEWIDTH, color=[1 2 3],
                   labels=reshape([STATELABELS...], 1, 6), ylabel="State\n(C/N mol)", 
                   xlabel="")
    rateplot = plot(ustrip.(hcat(model.records[1].vars.rate, model.records[2].vars.rate));
                    ylabel="Growth rate\n(mol mol^-1 d^-1)", xlabel="Month",
                    labels=["Shoot" "Root"], linewidth=LINEWIDTH, alpha=alpha)
    plt = plot(solplot, rateplot; xlims=(1,ustrip(tstop)), xticks=xticks, link=:x,
               layout=grid(2, 1, heights=[0.6, 0.4]), size=(600,400), dpi=DPI)
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
    xticks=(ustrip(0hr:MONTH_HOURS:tstop), MONTHS[1:months+1])
    n = length(u) ÷ 2 
    solt = sol'
    solt[:, 4:6] = solt[:, 4:6] * -1
    solplot = plot(solt, legend=:topleft, linewidth=LINEWIDTH, color=[1 2 3],
                   labels=reshape([STATELABELS...], 1, 6), ylabel="State\n(C/N mol)", 
                   xticks=xticks, xlabel="")
    assimplot = plot(ustrip.(hcat(model.records[1].J[2,1,:])); 
                     ylabel="Assimilation\n(C-mol hr^-1)", xticks=xticks, linewidth=LINEWIDTH, legend=false)
    swpplot = plot(hcat(model.records[1].vars.swp); yaxis=("Maximum\nsoil water\npotential"), 
                   xlabel="Month", xticks=xticks, linewidth=LINEWIDTH, 
                   labels=RANGE, legend=false)
    plt = plot(solplot, assimplot, swpplot; xlims=(1,ustrip(tstop)), link=:x,
         layout=grid(3, 1, heights=[0.5, 0.25, 0.25]), size=(600,600), dpi=DPI)
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
    xticks=(ustrip(0hr:MONTH_HOURS:tstop), MONTHS[1:months])
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
    plot(plt; xlims=(1, ustrip(tstop)), size=(600,600), dpi=DPI, legend=:none)
    savefig("plots/scaling")
    plt
end

function plot_tempcorr(model, environment, u, envstart, months)
    tstop = round(typeof(1hr), months * MONTH_HOURS)
    rnge = ustrip(envstart:1hr:envstart+tstop)
    plotstart = 1hr
    model.environment = environment
    model.environment_start[] = envstart
    model.dead[] = false
    tspan = ustrip.((1, tstop))
    prob = DiscreteProblem(model, ustrip(u), tspan)
    sol = solve(prob, FunctionMap(scale_by_time = true))
    n = length(u) ÷ 2 
    xticks=(ustrip(0hr:MONTH_HOURS:tstop), MONTHS[1:months+1])

    swp = -1 .* ustrip(soilwaterpotential(model.environment)) # .|> MPa
    st = soiltemperature(model.environment) .|> °C
    vpd = vapour_pressure_deficit.(airtemperature(model.environment), relhumidity(model.environment))
    rad = radiation(model.environment)
    soilcolors = permutedims(reverse(get.(Ref(ColorSchemes.copper), 0.0:1/7:1)))

    radplot = plot(rad[rnge]; ylab="Radiation", color=:black, legend=:none)
    soilwaterplot = plot(swp[rnge, :]; yflip=true, yscale=:log10, ylab=ENVLABS[1], color=soilcolors, 
                         labels=INCREMENTS, legend=:best, linewidth=LINEWIDTH, legendfontsize=7, alpha=alpha)
    soiltempplot = plot(st[rnge, :]; ylab=ENVLABS[3], ylims=(-10, 80), 
                        color=soilcolors, linewidth=LINEWIDTH, legend=false, alpha=alpha)
    vpdplot = plot(vpd[rnge, :]; ylab=ENVLABS[2], alpha=alpha, color=[:blue :green], 
                  labels=RANGE, linewidth=LINEWIDTH, legend=:bottom, legendfontsize=7)
    tempcorrplot = plot(ustrip.(hcat(model.records[1].vars.tempcorrection, model.records[2].vars.tempcorrection));
                    ylabel="Temperatuire\ncorrection", xlabel="", linewidth=LINEWIDTH, alpha=alpha, 
                    labels=["Shoot" "Root"], legend=:bottom, legendfontsize=7)
    rateplot = plot(ustrip.(hcat(model.records[1].vars.rate, model.records[2].vars.rate));
                    ylabel="Growth\nrate\n(d^-1)", xlabel="Month", 
                    labels=["Shoot" "Root"], linewidth=LINEWIDTH, alpha=alpha)

    plts = (radplot, soilwaterplot, soiltempplot, vpdplot, tempcorrplot, rateplot)
    plt = plot(plts...; xticks=xticks, xlims=(1, ustrip(tstop)), linewidth=LINEWIDTH, link=:x, layout=grid(length(plts), 1), 
               size=(600,600), dpi=DPI)
    savefig("plots/tempcorr")
    plt
end

theme(:sand)
theme(:solarized)
theme(:juno)
theme(:solarized_light)
theme(:wong2)
x = 1.5

pyplot()
# pyplot()
# plotly()

environments, _ = loadenvironments(dir)
models = OrderedDict()
modeldir = joinpath(dir, "models")
include.(joinpath.(Ref(modeldir), readdir(modeldir)));
envstart = round(typeof(1hr), 7 * MONTH_HOURS)
model = models[:bbiso]
# model = models[:bb]
environment = environments[:t1]
INCREMENTS = reshape([string.(Microclimate.get_increments(environment))...], 1, 8)
RANGE = reshape([string.(Microclimate.get_range(environment))...], 1, 2)
u = states["Large seed"]

crossover = plot_crossover(model, environment, u, envstart, 1);
crossover
assim = plot_assim(model, environment, u, envstart, 6);
tempcorr = plot_tempcorr(deepcopy(model), environment, u, envstart, 2);
tempcorr
scaling = plot_scaling(deepcopy(models[:bb]), environment, u, envstart, 12, 1.5);

gr()
envstart = 1hr
allplot = plot_years(model, environments, states, envstart)
