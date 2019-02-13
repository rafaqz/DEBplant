using Revise, Unitful, Flatten, FieldMetadata, OrdinaryDiffEq
using Photosynthesis, Microclimate, DynamicEnergyBudgets
using InteractBulma, InteractBase, Blink, WebIO, Observables, CSSUtil, Mux
using Plots, UnitfulPlots, PlotNested, StatsPlots

import Plots:px, pct, GridLayout
using Photosynthesis: potential_dependence
using DynamicEnergyBudgets: STATE, STATE1, TRANS, TRANS1, shape_correction, define_organs,
      photosynthesis, split_state, HasCN, HasCNE, has_reserves, tempcorr_pars,
      assimilation_pars, assimilation_vars, parconv, w_V
using Unitful: °C, K, Pa, kPa, MPa, J, kJ, W, L, g, kg, cm, m, s, hr, d, mol, mmol, μmol, σ

using JLD2


if "DEBSCRIPTS" in keys(ENV)
    dir = ENV["DEBSCRIPTS"]
else
    dir = pwd()
end

include(joinpath(dir, "util/hide.jl"))


# Import environment locations
locationspath = joinpath(dir, "microclimate/locations.jld")
@load locationspath tas desert qld

environment = tas
# environment = desert
# environment = qld

plotly()

const tspan = (0:1:length(environment.radiation) - 1) * hr
const statelabels = tuple(vcat([string("Shoot ", s) for s in STATE], [string("Root ", s) for s in STATE])...)
global uimodel = DynamicEnergyBudgets.PlantCN(environment=environment, time=tspan, environment_start=1u"hr")
global uisol

init_state(modelobs::Observable) = init_state(modelobs[])
init_state(model) = init_state(has_reserves.(define_organs(model, 1hr)))
init_state(::NTuple{2,HasCN}) = [0.0, 1e-4, 0.0, 1e-4, 1e-4, 1e-4, 0.0, 1e-4, 0.0, 0.01, 0.0005, 0.0]mol
init_state(::NTuple{2,HasCNE}) = [0.0, 1e-4, 0.0, 1e-4, 1e-4, 1e-4, 0.0, 1e-4, 0.0, 1e-4, 1e-4, 0.01]mol
init_state(::NTuple{3,HasCN}) = [0.0, 1e-4, 0.0, 1e-4, 1e-4, 1e-4, 0.0, 1e-4, 0.0, 1e-4, 1e-4, 1e-4, 0.0, 1e-4, 0.0, 0.01, 0.002, 0.0]mol
init_state(::NTuple{3,HasCNE}) = [0.0, 1e-4, 0.0, 1e-4, 1e-4, 1e-4, 0.0, 1e-4, 0.0, 1e-4, 1e-4, 1e-4, 0.0, 1e-4, 0.0, 1e-4, 1e-4, 0.01]mol

photo(o, u, x) = begin
    update_photovars(assimilation_vars(o), x)
    photosynthesis(assimilation_pars(o), o, u) / (w_V(o) * assimilation_pars(o).SLA)
end

pot(o, x) = potential_dependence(assimilation_pars(o).potential_modifier, x)

update_photovars(v::CarbonVars, x) = begin
    v.J_L_F = x*parconv
    v.soilwaterpotential = zero(v.soilwaterpotential)
end
update_photovars(v::Photosynthesis.PhotoVars, x) = v.par = x * oneunit(v.par)

function sol_plot(model::AbstractOrganism, params::AbstractVector, u::AbstractVector, tstop::Number, envstart::Number)
    length(params) > 0 || return (plot(), plot())
    m2 = reconstruct(model, params)

    model.params = m2.params
    model.shared = m2.shared
    model.records = m2.records
    model.environment_start[] = envstart
    model.dead[] = false

    global uimodel = model

    println("u, tstop: ", (u, tstop))
    prob = DiscreteProblem(model, u, (oneunit(tstop), tstop))
    sol = solve(prob, FunctionMap(scale_by_time = true))
    background_color = model.dead[] ? :red : :white
    # solplot1 = plot(sol, vars = [1:6...], plotdensity=400, legend=:topleft, background_color=background_color,
    #                 labels=reshape(statelabels[1:6], 1, 6), ylabel="State (CMol)",
    #                 xlabel=string(model.params[1].name, " : ",  typeof(model.params[1].assimilation_pars).name, " - time (hr)"))
    # solplot2 = plot(sol, vars = [7:12...], plotdensity=400, legend=:topleft, background_color=background_color,
    #                 labels=reshape(statelabels[7:12], 1, 6), ylabel="State (CMol)",
    #                 xlabel=string(model.params[2].name, " : ", typeof(model.params[2].assimilation_pars).name, " - time (hr)"))
    # plot(solplot1, solplot2, layout=Plots.GridLayout(2, 1))
    s = sol' .* model.shared.core_pars.w_V
    s1 = s[:, 1:6]
    s2 = s[:, 7:12]
    global uisol = s
    solplot = plot(sol.t, s1, labels=reshape([statelabels[1:6]...], 1, 6)), plot(sol.t, s2, labels=reshape([statelabels[7:12]...], 1, 6))
    # if plotarea
    # organs = define_organs(model, 1hr)
    # o = organs[1]
        # areaplot = plot(s[:, 2] .* assimilation_pars(o).SLA, ylabel="C uptake modification", xlabel="Surface area")
        # plot(solplot, areaplot, layout=Plots.GridLayout(2, 1))
    # else
        # solplot
    # end
end


function make_plot(model::AbstractOrganism, params::AbstractVector, u::AbstractVector, solplots, vars, env, flux,
                   plottemp::Bool, plotscale::Bool, plotphoto::Bool, plotpot::Bool, tstop::Number)
    plotsize=(1700, 800)

    if length(params) > 0
        model = reconstruct(model, params)
    end
    envstart = model.environment_start[]
    tspan = 1:ustrip(tstop)
    varplots = plot_selected(model.records, vars, tspan)
    envplots = plot_selected(model.environment, env, ustrip(envstart):ustrip(envstart+tstop))
    fluxplots = []
    if length(flux) > 0
        plot_fluxes!(fluxplots, flux[1], model.records[1].J, tspan)
        plot_fluxes!(fluxplots, flux[2], model.records[1].J1, tspan)
        plot_fluxes!(fluxplots, flux[3], model.records[2].J, tspan)
        plot_fluxes!(fluxplots, flux[4], model.records[2].J1, tspan)
    end
    timeplots = plot(solplots..., varplots..., envplots..., fluxplots...,
                     layout=Plots.GridLayout(length(solplots)+length(varplots)+length(envplots)+length(fluxplots), 1))

    organs = define_organs(model, 1hr)
    o = organs[1]
    state = split_state(organs, u, 0)

    subplots = []

    plottemp && push!(subplots, plot(x -> tempcorr(x, tempcorr_pars(o)), 0.0°C:1.0°C:50.0°C,
                                     legend=false, ylabel="Correction", xlabel="°C"))
    plotscale && push!.(Ref(subplots), plot_shape.(model.params))
    plotphoto && push!(subplots, plot(x -> photo(o, state[1], x), (0*W*m^-2:10*W*m^-2:1000*W*m^-2),
                                      ylabel="C uptake", xlabel="Irradiance"))
    plotpot && push!(subplots, plot(x -> pot(o, x), (0.0kPa:-1.0kPa:-5000kPa),
                                    ylabel="C uptake modification", xlabel="Soil water potential"))

    if length(subplots) > 0
        funcplots = plot(subplots..., layout=Plots.GridLayout(length(subplots), 1))
        l = Plots.GridLayout(1, 2)
        l[1, 1] = GridLayout(1, 1, width=0.8pct)
        l[1, 2] = GridLayout(1, 1, width=0.2pct)
        plot(timeplots, funcplots, size=plotsize, layout=l)
    else
        plot(timeplots, size=plotsize)
    end
end

plot_shape(p) = plot(x -> shape_correction(p.shape_pars, x), (0.0mol:0.01mol:100.0mol),
                       legend=false, ylabel="Correction", xlabel="CMols Stucture")


plot_fluxes!(plots, obs, J, tspan) = begin
    ps = []
    labels = []
    for y in 1:size(J, 1), x in 1:size(J, 2)
        if obs[y, x]
            push!(labels, "$y, $x")
            push!(ps, map(t -> ustrip(J[y, x, t]), tspan .+ oneunit(eltype(tspan))))
        end
    end
    if length(ps) > 0
        push!(plots, plot(ps, labels=reshape(labels, 1, length(labels))))
    end
end

function spreadwidgets(widgets; cols = 5)
    vboxes = []
    widget_col = []
    colsize = ceil(Int, length(widgets)/cols)
    # Build vbox columns for widgets
    for i = 1:length(widgets)
        push!(widget_col, widgets[i])
        if rem(i, colsize) == 0
            push!(vboxes, vbox(widget_col...))
            widget_col = []
        end
    end
    # Push remaining widgets
    push!(vboxes, vbox(widget_col...))
    hbox(vboxes...)
end

allsubtypes(x) = begin
    st = subtypes(x)
    if length(st) > 0
        allsubtypes((st...,))
    else
        (x,)
    end
end
allsubtypes(t::Tuple{X,Vararg}) where X = (allsubtypes(t[1])..., allsubtypes(Base.tail(t))...,)
allsubtypes(::Tuple{}) = ()

assimvars(::AbstractCAssim) = DynamicEnergyBudgets.ShootVars()
assimvars(::AbstractNAssim) = DynamicEnergyBudgets.RootVars()
assimvars(::FvCBPhotosynthesis) = DynamicEnergyBudgets.FvCBShootVars()

build_model(params, su, cat, rate, prod, mat, rej, trans, shape, allometry, assim1, assim2,
            temp, feedback, env, envstart) = begin
    pkwargs = (maturity_pars=mat(), rate_formula=rate(),
               production_pars=prod(), rejection_pars=rej(), trans_pars=trans(),
               shape_pars=shape(), allometry_pars=allometry())
    p1 = params(;pkwargs..., assimilation_pars=assim1())
    p2 = params(;deepcopy(pkwargs)..., assimilation_pars=assim2())
    sh = SharedParams(su_pars=su(), tempcorr_pars=temp(), catabolism_pars=cat(), feedback_pars=feedback())
    envstart = setindex!(Array{typeof(envstart),0}(undef), envstart)
    Plant(params=(p1, p2), shared=sh, environment=environment, time=tspan, environment_start=envstart,
             vars=(assimvars(assim1()), assimvars(assim2())))
end

function muxapp(req) # an "App" takes a request, returns the output

    envobs = Observable{Any}(environment);

    tstoptext = textbox(value="8760", label="Timespan")
    envstart = slider(1hr:1hr:tspan.stop, value = 1hr, label="Environment start time")
    envdrop = dropdown([:tas, :desert, :qld], value=:tas, label="Environment")

    plottemp = checkbox("Plot TempCorr")
    plotscale = checkbox("Plot Scaling Function")
    plotphoto = checkbox("Plot Photosynthesis")
    plotpot = checkbox("Plot Potential dependence")

    controlbox = hbox(tstoptext, envstart, envdrop, plottemp, plotscale, plotphoto, plotpot)

    reload = button("Reload")
    paramsdrop = dropdown([allsubtypes(AbstractParams)...], value=Params, label="Params")
    sudrop = dropdown([allsubtypes(AbstractSU)...], value=ParallelComplementarySU, label="SU")
    ratedrop = dropdown([allsubtypes(AbstractRate)...], label="Rate")
    catabolismdrop = dropdown([allsubtypes(AbstractCatabolism)...], label="Catabolism")
    productiondrop = dropdown([Nothing, allsubtypes(AbstractProduction)...], label="Production")
    maturitydrop = dropdown([Nothing, allsubtypes(AbstractMaturity)...], label="Maturity")
    rejectiondrop = dropdown([Nothing, allsubtypes(AbstractRejection)...], value=LosslessRejection, label="Rejection")
    translocationdrop = dropdown([Nothing, allsubtypes(AbstractTranslocation)...], label="Translocation")
    shapedrop = dropdown([allsubtypes(AbstractShape)...], value=Isomorph, label="Shape")
    allometrydrop = dropdown([Nothing, allsubtypes(AbstractAllometry)...], value=SqrtAllometry, label="Allometry")
    assimdrop1 = dropdown([Nothing, allsubtypes(AbstractAssim)...], value=KooijmanWaterPotentialPhotosynthesis, label="Assimilation")
    assimdrop2 = dropdown([Nothing, allsubtypes(AbstractAssim)...], value=ConstantNAssim, label="Assimilation")
    tempdrop = dropdown([Nothing, allsubtypes(AbstractTempCorr)...], value=TempCorrLowerUpper, label="Temp Correction")
    feedbackdrop = dropdown([Nothing, allsubtypes(AbstractStateFeedback)...], value=DissipativeAutophagy, label="State Feedback")
    dropbox = vbox(hbox(reload, paramsdrop, sudrop, catabolismdrop, ratedrop, productiondrop, maturitydrop, rejectiondrop, translocationdrop),
                   hbox(shapedrop, allometrydrop, assimdrop1, assimdrop2, tempdrop, feedbackdrop))

    solplotobs = Observable{Any}([plot(), plot()])
    plotobs = Observable{Any}(plot())



    modelobs = Observable{Any}(DynamicEnergyBudgets.PlantCN(environment=environment, time=tspan, environment_start=1u"hr"))
    map!(modelobs, envdrop) do e
        if e == :tas
            env = tas
        elseif e == :desert
            env = desert
        elseif e == :qld
            env = qld
        end
        modelobs[].environment = env
        modelobs[]
    end

    reload_model() = begin
        modelobs[] = build_model(paramsdrop[], sudrop[], catabolismdrop[], ratedrop[], productiondrop[], maturitydrop[], rejectiondrop[], translocationdrop[], shapedrop[],
                                 allometrydrop[], assimdrop1[], assimdrop2[], tempdrop[], feedbackdrop[], environment, envstart[])
        global uimodel = modelobs[]
        # plotobs[] = make_plot(modelobs[], varobs[], stateobs[], paramobs[], fluxobs[], observe(tstop)[])
    end

    on(observe(reload)) do x
        reload_model()
    end


    paramsliders = Observable{Vector{Widget{:slider}}}(Widget{:slider}[]);
    paramobs = Observable{Vector{Any}}([]);
    parambox = Observable{typeof(dom"div"())}(dom"div"());

    varchecks = Observable{Vector{Widget{:checkbox}}}(Widget{:checkbox}[]);
    varobs = Observable{Vector{Bool}}(Bool[]);
    varbox = Observable{typeof(dom"div"())}(dom"div"());

    envchecks = Observable{Vector{Widget{:checkbox}}}(Widget{:checkbox}[]);
    envobs = Observable{Vector{Bool}}(Bool[]);
    envbox = Observable{typeof(dom"div"())}(dom"div"());

    state = init_state(modelobs[])
    statesliders = Observable{Vector{Widget{:slider}}}(Widget{:slider}[]);
    stateobs = Observable{typeof(state)}(state);
    statebox = Observable{typeof(dom"div"())}(dom"div"());

    make_varchecks(model) = begin
        checks = plotchecks(model.records)
        map!((x...) -> [x...], varobs, observe.(checks)...)
        checks
    end

    make_envchecks(model) = begin
        checks = plotchecks(model.environment)
        map!((x...) -> [x...], envobs, observe.(checks)...)
        checks
    end

    make_paramsliders(m) = begin
        params = flatten(Vector, m)
        fnames = fieldnameflatten(Vector, m)
        parents = parentflatten(Vector, m)
        limits = metaflatten(Vector, m, FieldMetadata.limits)
        descriptions = metaflatten(Vector, m, FieldMetadata.description)
        units = metaflatten(Vector, m, FieldMetadata.units)
        log = metaflatten(Vector, m, FieldMetadata.logscaled)
        attributes = broadcast((p, n, d, u) -> Dict(:title => "$p.$n: $d $u"), parents, fnames, descriptions, units)
        sl = broadcast(limits, fnames, params, attributes, log) do l, n, p, a, lg
            println((l, n, p, a, lg))
            # Use a log range if specified in metadata
            rnge = lg ? vcat(exp10.(range(log10(abs(ustrip(l[1]))), stop=sign(log10(abs(ustrip(l[2])))) , length=100) * sign(l[2])) * l[2]) : collect(l[1]:(l[2]-l[1])/100:l[2])
            # rnge = l[1]:(l[2]-l[1])/100:l[2]
            InteractBase.slider(rnge, label=string(n), value=p, attributes=a)
        end
        map!((x...) -> [x...], paramobs, throttle.(0.2, observe.(sl))...)
        sl
    end

    make_statesliders(m) = begin
        sl = state_slider.(statelabels, init_state(m))
        map!((x...) -> [x...], stateobs, throttle.(0.2, observe.(sl))...)
        sl
    end

    flux_grids = (make_grid(STATE, TRANS), make_grid(STATE1, TRANS1,),
                  make_grid(STATE, TRANS), make_grid(STATE1, TRANS1));
    fluxbox = hbox(arrange_grid.(flux_grids)...)
    fluxobs = map((g...) -> g, observe_grid.(flux_grids)...)
    tstopobs = map(throttle(0.2, observe(tstoptext))) do t
        int = tryparse(Int, t)
        int == nothing ? 8760hr : int * hr
    end
    envstartobs = throttle(0.2, observe(envstart))

    map!(make_varchecks, varchecks, modelobs)
    map!(make_envchecks, envchecks, modelobs)
    map!(make_statesliders, statesliders, modelobs)
    map!(make_paramsliders, paramsliders, modelobs)
    map!(sol_plot, solplotobs, modelobs, paramobs, stateobs, tstopobs, envstartobs)
    map!(make_plot, plotobs, modelobs, paramobs, stateobs, solplotobs, varobs, envobs, fluxobs, observe(plottemp), observe(plotscale), observe(plotphoto), observe(plotpot), tstopobs)

    map!(c -> vbox(hbox(c...)), varbox, varchecks)
    halfenv = min(12, length(envchecks[]))
    map!(c -> vbox(hbox(c[1:halfenv]...), hbox(c[halfenv+1:end]...)), envbox, envchecks)
    map!(s -> spreadwidgets(s), parambox, paramsliders)
    map!(s -> spreadwidgets(s), statebox, statesliders)

    reload_model()

    ui = vbox(controlbox, dropbox, fluxbox, plotobs, varbox, envbox, parambox, statebox);
    # dom"div"()
end

state_slider(label, val) =
    slider(vcat(exp10.(range(-4, stop = 1, length = 100))) * mol, label=label, value=val)

make_grid(rownames, colnames) = begin
    rows = length(rownames)
    cols = length(colnames)
    [checkbox(false, label = join([string(rownames[r]), string(colnames[c])], ",")) for r = 1:rows, c = 1:cols]
end

arrange_grid(a) = hbox((vbox(a[:,col]...) for col in 1:size(a, 2))...)

observe_grid(a) = map((t...) -> [t[i + size(a,1) * (j - 1)] for i in 1:size(a,1), j in 1:size(a,2)], observe.(a)...)

ui = muxapp(nothing)
w = Window(Dict("webPreferences"=>Dict("zoomFactor"=>0.6)));
Blink.AtomShell.@dot w webContents.setZoomFactor(0.6)
body!(w, ui);
# opentools(w)

# webio_serve(page("/", req -> muxapp(req)), 8000)
# @relimits @redefault struct Maintenance
    # j_E_mai   | 0.01 | [1e-4, 1.0]
# end
