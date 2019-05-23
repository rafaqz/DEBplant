dir = "DEBSCRIPTS" in keys(ENV) ? ENV["DEBSCRIPTS"] : pwd()
include(joinpath(dir, "load.jl"))
include(joinpath(dir, "plantstates.jl"))

using PlotThemes, LaTeXStrings

MONTH_HOURS = 365.25 / 12 * 24hr
YLIMS = (-3,7)
LIFESPAN = 4380
LINEWIDTH = 1.6
MULTIYEARLINEWIDTH = 0.2
ENVLABS = ("Soil water\npotential", "Relative\nhumidity", "Soil\ntemperature")
DPI = 130

import Base: round
round(::Type{T}, x::Quantity) where {T<:Quantity} = T(round(typeof(one(T)), uconvert(unit(T), x).val)) 

function plot_microclim(model, envname, rnge = 1:1:size(radiation(model.environment), 1); 
                        ylabs=ENVLABS, kwargs...)
    swp = -1 .* ustrip(soilwaterpotential(model.environment)) # .|> MPa
    st = soiltemperature(model.environment) .|> °C
    rh = relhumidity(model.environment)
    soilcolors = permutedims(reverse(get.(Ref(ColorSchemes.copper), 0.0:1/7:1)))
    # yticks = ["-10⁰ KPa", "-10¹ KPa", "-10² KPa", "-10³ KPa", "-10⁴ KPa"]
    # yticks = [-10^0u"kPa", -10^1u"kPa", -10^2u"kPa", -10^3u"kPa", -10^4u"kPa"]
    yticks = ([10^0, 10^1, 10^2, 10^3, 10^4], 
              string.(Ref("-10"), [L"_0", L"_1", L"_2", L"_3", L"_4"], Ref(" kPa")))
    soilwaterplot = plot(rnge * hr, swp[rnge, :]; yflip=true, yscale=:log10, ylab=ylabs[1], 
                         color=soilcolors, legend=:none, xshowaxis=false, yticks=yticks, kwargs...)
    rhplot = plot(rnge * hr, rh[rnge, :]; ylab=ylabs[2], linealpha=0.8,
                  legend=:none, color=[:black :darkgrey], xshowaxis=false, kwargs...)
    soiltempplot = plot(rnge * hr, st[rnge, :]; xlab=titlecase(string(envname)), ylab=ylabs[3],
                        ylims=(-10, 80), color=soilcolors, legend=:none, kwargs...)
    soilwaterplot, rhplot, soiltempplot
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
    plot!(plt, rng, shootvals, linecolor=:black, xshowaxis=false)
    plot!(plt, rng, rootvals, linecolor=:grey, xshowaxis=false)
    return model, sol
end

" Plot growth starting every moth in the available timespan"
function plot_months(title, model, u, envstart; kwargs...)
    solplots = plot(; ylab=title, ylims=YLIMS,
                    legend=:none, link=:x, kwargs...)
    for i in 1:119
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
                push!(stateplots, plot_months("Shoot and Root\nmass", model, states[statename],  envstart))
            else
                push!(stateplots, plot_months("", model, states[statename], envstart; yshowaxis=false))
            end
        end
        # xaxis = ((0,8760*11), 0:8760*2:8760*11)
        xticks = (collect(0:8760*2:8760*11), string.(collect(2001:2:2011)))
        if i == 1
            microplots = plot_microclim(model, envname; linewidth=MULTIYEARLINEWIDTH)
            heights, num_plots = calc_heights(stateplots, microplots)
            push!(locplots, plot(stateplots..., microplots...;
                 xticks=xticks, link=:x, margin=0px, right_margin=0px,
                 layout=grid(num_plots, 1, heights=heights)))
        else
            microplots = plot_microclim(model, envname; linewidth=MULTIYEARLINEWIDTH, ylabs=("", "", ""), yshowaxis=false)
            heights, num_plots = calc_heights(stateplots, microplots)
            push!(locplots, plot(stateplots..., microplots...;
                 xticks=xticks, link=:x, margin=0px, left_margin=-50px,
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
    model.environment = environments[:t1]
    solplot = plot(yaxis=(name), legend=:none, link=:x)
    model, sol = plot_sol!(solplot, model, states[name], envstart)
    rateplot = plot(hcat(model.records[1].vars.rate, model.records[2].vars.rate); 
                    yaxis=("Growth rate"), labels=["Shoot" "Root"])
    assimplot = plot(model.records[1].J[4,1,:]; yaxis=("Assimilation"))
    plt = plot(solplot, rateplot, assimplot, microplots...; xaxis=((0,LIFESPAN), 0:500:4000), 
               link=:x, layout=grid(6, 1, heights=[0.3, 0.2, 0.15, 0.15, 0.1, 0.1]), size=(600,600), dpi=DPI)
    savefig("plots/single")
    plt
end

function plot_crossover(model, environment, u, envstart, tstop)
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
    # solplot1 = plot(sol, tspan=tspan, vars = [1:n...], plotdensity=400, legend=:topleft,
    #                 labels=reshape([STATELABELS[1:n]...], 1, n), 
    #                 ylabel="Shoot state\n(C/N-mol)", xlabel="", linewidth=1.5)
    # solplot2 = plot(sol, tspan=tspan, vars = [n+1:2n...], plotdensity=400, legend=:topleft,
    #                 labels=reshape([STATELABELS[n+1:2n]...], 1, n), 
    #                 ylabel="Root state\n(C/N-mol)", xlabel="", linewidth=1.5)
    solplot = plot(solt, legend=:topleft, linewidth=LINEWIDTH, color=[1 2 3],
                   labels=reshape([STATELABELS...], 1, 6), ylabel="State\n(C/N mol)", xlabel="")
    rateplot = plot(ustrip.(hcat(model.records[1].vars.rate, model.records[2].vars.rate));
                    ylabel="Growth rate\n(mol mol^-1 d^-1)", xlabel="Time (hr)", 
                    labels=["Shoot" "Root"], linewidth=LINEWIDTH)
    plt = plot(solplot, rateplot; xlims=(1,tstop), link=:x,
         layout=grid(2, 1, heights=[0.6, 0.4]), size=(600,400), dpi=DPI)
    savefig("plots/crossover")
    plt
end

function plot_assim(model, environment, u, envstart, tstop)
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
    solplot = plot(solt, legend=:topleft, linewidth=LINEWIDTH, color=[1 2 3],
                   labels=reshape([STATELABELS...], 1, 6), ylabel="State\n(C/N mol)", xlabel="")
    # solplot2 = plot(solt, vars = [n+1:2n...], legend=:topleft, linewidth=LINEWIDTH,
                    # labels=reshape([STATELABELS[n+1:2n]...], 1, n), ylabel="Root state\n(C/N mol)", xlabel="")
    assimplot = plot(ustrip.(hcat(model.records[1].J[2,1,:])); ylabel="Assimilation\n(C-mol hr^-1)", linewidth=LINEWIDTH, legend=false)
    swpplot = plot(hcat(model.records[1].vars.swp); yaxis=("Maximum\nsoil water\npotential"), 
                   xlabel="Time (hr)", linewidth=LINEWIDTH, legend=false)
    plt = plot(solplot, assimplot, swpplot; xlims=(1,tstop), link=:x,
         layout=grid(3, 1, heights=[0.5, 0.25, 0.25]), size=(600,600), dpi=DPI)
    savefig("plots/assim")
    plt
end

function plot_scaling(model, environment, u, envstart, tstop, x)
    plotstart = 1hr
    model.environment = environment
    model.environment_start[] = envstart
    model.dead[] = false
    tspan = ustrip.((plotstart, tstop))
    scaling = oneunit(model.params[2].shape_pars.M_Vscaling)
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
        annotate!(plt, tstop, solt[end, 1] - 0.05, text(string(round(ustrip(s), digits=2)), 7))
    end
    plot(plt; xlims=(1,tstop), size=(600,600), dpi=DPI, legend=:none)
    savefig("plots/scaling")
    plt
end

function plot_tempcorr(model, environment, u, envstart, tstop)
    rnge = ustrip(envstart):1:ustrip(envstart)+tstop
    plotstart = 1hr
    model.environment = environment
    model.environment_start[] = envstart
    model.dead[] = false
    tspan = ustrip.((1, tstop))
    prob = DiscreteProblem(model, ustrip(u), tspan)
    sol = solve(prob, FunctionMap(scale_by_time = true))
    n = length(u) ÷ 2 

    swp = -1 .* ustrip(soilwaterpotential(model.environment)) # .|> MPa
    st = soiltemperature(model.environment) .|> °C
    rh = relhumidity(model.environment)
    rad = radiation(model.environment)
    soilcolors = permutedims(reverse(get.(Ref(ColorSchemes.copper), 0.0:1/7:1)))

    radplot = plot(rad[rnge]; ylab="Radiation", color=:black, legend=:none)
    soilwaterplot = plot(swp[rnge, :]; yflip=true, yscale=:log10, ylab=ENVLABS[1], 
                         color=soilcolors, labels=string.(mclm_increments))
    rhplot = plot(rh[rnge, :]; labels=string.(mclm_range), ylab=ENVLABS[2], linealpha=0.8, 
                  color=[:blue :green], legend=:bottomright)
    soiltempplot = plot(st[rnge, :]; ylab=ENVLABS[3], ylims=(-10, 80), 
                        color=soilcolors, labels=string.(mclm_increments))
    tempcorrplot = plot(ustrip.(hcat(model.records[1].vars.tempcorrection, model.records[2].vars.tempcorrection));
                    ylabel="Temperatuire\ncorrection", xlabel="", 
                    labels=["Shoot" "Root"], linewidth=LINEWIDTH, legend=:bottomright)
    rateplot = plot(ustrip.(hcat(model.records[1].vars.rate, model.records[2].vars.rate));
                    ylabel="Growth\nrate\n(d^-1)", xlabel="Time (hr)", 
                    labels=["Shoot" "Root"], linewidth=LINEWIDTH)

    plts = (radplot, soilwaterplot, soiltempplot, rhplot, tempcorrplot, rateplot)
    plt = plot(plts...; xlims=(1,tstop), link=:x, layout=grid(length(plts), 1), size=(600,600), dpi=DPI)
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
envstart = 20000.0hr
model = models[:bbiso]
# model = models[:bb]
environment = environments[:t1]
u = states["Large seed"]
tstop = 1000
crossover = plot_crossover(model, environment, u, envstart, tstop);

tstop = 4000
assim = plot_assim(model, environment, u, envstart, tstop);

tstop = 1500
envstart = 1000hr
tempcorr = plot_tempcorr(deepcopy(model), environment, u, envstart, tstop);

tstop = 8760
envstart = 1hr
scaling = plot_scaling(deepcopy(models[:bb]), environment, u, envstart, tstop, 1.5);
tempcorr
mclm_increments = Microclimate.get_increments(environment)
mclm_range = Microclimate.get_range(environment)


# plt = plot_single(model, environments, states, 1000.0hr);
gr()
allplot = plot_years(model, environments, states, envstart)

# plt = plot([1,2,3]hr, 1:3)
# plot!(plt, [1,2,3]hr, 4:6)
# plot(plt, plot(rand(6)), link=:x, layout=(2,1))
