using Revise, Setfield, InteractBulma, InteractBase, Blink, WebIO, Observables, CSSUtil, Mux
using Plots, UnitfulPlots, StatsPlots, PlotNested, Codify, Select, ColorSchemes

include(joinpath(dir, "load.jl"))
include(joinpath(dir, "util/hide.jl"))

const STATELABELS = tuple(vcat([string("Shoot ", s) for s in STATE], [string("Root ", s) for s in STATE])...)

mutable struct ModelApp{M,E,T}
    models::M
    environments::E
    tspan::T
    savedmodel
end

init_state(modelobs::Observable) = init_state(modelobs[])
init_state(model) = init_state(has_reserves.(define_organs(model, 1hr)))
init_state(::NTuple{2,HasCN}) = [0.0, 1e-4, 0.0, 1e-4, 1e-4, 1e-4, 0.0, 1e-4, 0.0, 0.01, 0.0005, 0.0]mol
init_state(::NTuple{2,HasCNE}) = [0.0, 1e-4, 0.0, 1e-4, 1e-4, 1e-4, 0.0, 1e-4, 0.0, 1e-4, 1e-4, 0.01]mol

plotly()


pot(p, x) = potential_dependence(p.potential_modifier, x)

photo(o, u, x) = photo(assimilation_vars(o), o, u, x)
photo(v::CarbonVars, o, u, x) = begin
    v.J_L_F = x*parconv
    v.soilwaterpotential = zero(v.soilwaterpotential)
    photosynthesis(assimilation_pars(o), o, u) / (w_V(o) * assimilation_pars(o).SLA)
end
photo(v::Photosynthesis.PhotoVars, o, u, x) = begin
    v.rnet = x
    v.par = x * parconv
    run_enbal!(assimilation_pars(o).photoparams, v)
    v.aleaf * assimilation_pars(o).SLA * w_V(o)
end

function sol_plot(model::AbstractOrganism, params::AbstractVector, u::AbstractVector, tstop::Number, envstart::Number, app)

    length(params) > 0 || return (plot(), plot())
    m2 = ulreconstruct(model, params)

    model.params = m2.params
    model.shared = m2.shared
    model.records = m2.records
    model.environment_start[] = envstart
    model.dead[] = false
    set_allometry(model, u)
    # println("vars: ", model.records[1].vars.rate)

    app.savedmodel = model

    # println("u, tstop: ", (ustrip(u), ustrip(tstop)))
    # println("num params: ", length(params))
    prob = DiscreteProblem{true}(model, ustrip(u), (one(tstop), ustrip(tstop)))
    local sol
    try
        sol = solve(prob, FunctionMap(scale_by_time = true))
    catch e
        println(e)
        return plot(), plot()
    end
    solplot1 = plot(sol, vars = [1:6...], plotdensity=400, legend=:topleft,
                    labels=reshape([STATELABELS[1:6]...], 1, 6), ylabel="State (CMol)",
                    xlabel=string(model.params[1].name, " : ",  typeof(model.params[1].assimilation_pars).name, " - time (hr)"))
    solplot2 = plot(sol, vars = [7:12...], plotdensity=400, legend=:topleft,
                    labels=reshape([STATELABELS[7:12]...], 1, 6), ylabel="State (CMol)",
                    xlabel=string(model.params[2].name, " : ", typeof(model.params[2].assimilation_pars).name, " - time (hr)"))
    # plot(solplot1, solplot2, layout=Plots.GridLayout(2, 1))
    # s = sol' # .* model.shared.core_pars.w_V
    # s1 = view(s, :, 1:6)
    # s2 = s[:, 7:12]
    # solplot = (plot(sol.t, s1, labels=reshape([STATELABELS[1:6]...], 1, 6)),)
    # plot(sol.t, s2, labels=reshape([STATELABELS[7:12]...], 1, 6))
    # if plotarea
    # organs = define_organs(model, 1hr)
    # o = organs[1]
        # areaplot = plot(s[:, 2] .* assimilation_pars(o).SLA, ylabel="C uptake modification", xlabel="Surface area")
        # plot(solplot, areaplot, layout=Plots.GridLayout(2, 1))
    # else
        # solplot
    # end
    solplot1, solplot2
end


function make_plot(model::AbstractOrganism, u::AbstractVector, solplots, vars, env, flux,
                   plottemp::Bool, plotscale::Bool, plotphoto::Bool, plotpot::Bool, tstop::Number)
    plotsize=(1700, 800)

    envstart = model.environment_start[]
    tspan = 2:ustrip(tstop)
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

    plottemp && push!(subplots, plot(x -> tempcorr(x, tempcorr_pars(o)), K(0.0°C):1.0K:K(50.0°C),
                                     legend=false, ylabel="Correction", xlabel="°C"))
    plotscale && push!.(Ref(subplots), plot_shape.(model.params))
    plotphoto && push!(subplots, plot(x -> photo(o, state[1], x), (0*W*m^-2:5*W*m^-2:1000*W*m^-2),
                                      ylabel="C uptake", xlabel="Irradiance"))
    plotpot && push!(subplots, plot(x -> pot(assimilation_pars(o), x), (0.0kPa:-10.0kPa:-5000kPa),
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
allsubtypes(::Nothing) = (Nothing,)
allsubtypes(t::Tuple{X,Vararg}) where X = (allsubtypes(t[1])..., allsubtypes(Base.tail(t))...,)
allsubtypes(::Tuple{}) = ()
allsubtypes(t::Union) = (allsubtypes(t.a)..., allsubtypes(t.b)...,)

runplotly() = plotly()
rungr() = gr()

gr()


function (app::ModelApp)(req) # an "App" takes a request, returns the output

    throt = 0.1
    env = first(values(app.environments))
    envobs = Observable{Any}(env);


    tstoptext = textbox(value="8760", label="Timespan")
    envstart = slider(1hr:1hr:tspan.stop, value = 1hr, label="Environment start time")
    envdrop = dropdown(app.environments, value=env, label="Environment")
    modeldrop = dropdown(app.models, value=app.models[:init], label="Model")
    modelobs = Observable{Any}(first(values(app.models)))
    gr = button("GR")
    plotly = button("Plotly")
    save = button("Save")
    savename = textbox(value="modelname", label="Save name")

    plottemp = checkbox("Plot TempCorr")
    plotscale = checkbox("Plot Scaling Function")
    plotphoto = checkbox("Plot Photosynthesis")
    plotpot = checkbox("Plot Potential dependence")

    controlbox = hbox(save, savename, modeldrop, envdrop, tstoptext, envstart, plottemp, plotscale, plotphoto, plotpot, plotly, gr)

    reload = button("Reload")

    on(p -> rungr(), observe(gr))
    on(p -> runplotly(), observe(plotly))


    selectables = select(modelobs[])
    drops = []
    for (typ, def, lab) in selectables
        push!(drops, dropdown([allsubtypes(typ)...], value=def, label=lab))
    end
    half = length(drops) ÷ 2
    dropbox = vbox(hbox(drops[1:half]...), hbox(drops[half+1:end]...))
    getindex.(drops)

    solplotobs = Observable{Any}([plot(), plot()])
    plotobs = Observable{Any}(plot())

    map!(modelobs, observe(modeldrop), observe(envdrop)) do m, e
        load_model(m, e)
    end

    reload_model(app, drops) = begin
        m = updateselected(modelobs[], getindex.(drops))
        modelobs[] = load_model(m, envdrop[])
    end

    load_model(m, e) = begin
        m = update_vars(m, 87600)
        m.environment = e
        # Update all the component dropdowns
        setindex!.(drops, getindex.(select(m), 2))
        app.savedmodel = m
        m
    end


    on(observe(reload)) do x
        reload_model(app, drops)
    end

    on(observe(save)) do x
        savecode(app, savename[])
    end

    paramsliders = Observable{Vector{Any}}(Widget{:slider}[]);
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

    make_paramsliders(model) = begin
        params = ulflatten(Vector, model)
        fnames = fieldnameflatten(Vector, model)
        parents = parentflatten(Vector, model)
        backcolors = hex.(colorforeach(parents))
        limits = metaflatten(Vector, model, FieldMetadata.limits)
        descriptions = metaflatten(Vector, model, FieldMetadata.description)
        unts = metaflatten(Vector, model, FieldMetadata.units)
        log = metaflatten(Vector, model, FieldMetadata.logscaled)
        attributes = broadcast((p, b, n, d, u) -> Dict(:title => "$p.$n: $d $(u == nothing ? "" : u)"), parents, backcolors, fnames, descriptions, unts)
        sl = broadcast(limits, fnames, params, attributes, log, unts) do l, n, p, a, lg, u
            # println((l, n, p, a, lg))
            # Use a log range if specified in metadata
            rnge = lg ? vcat(exp10.(range(log10(abs(l[1])), stop=sign(log10(abs(l[2]))) , length=100) * sign(l[2])) * l[2]) : collect(l[1]:(l[2]-l[1])/100:l[2])
            InteractBase.slider(rnge, label=string(n), value=p, attributes=a)
        end
        paramobs[] = [s[] for s in sl]
        map!((x...) -> [x...], paramobs, throttle.(throt, observe.(sl))...)

        solplotobs[] = sol_plot(modelobs[], paramobs[], stateobs[], tstopobs[], envstartobs[], app)
        [dom"div[style=background-color:#$(backcolors[i])]"(sl) for (i, sl) in enumerate(sl)]
    end

    make_statesliders(m) = begin
        sl = state_slider.(STATELABELS, init_state(m))
        map!((x...) -> [x...], stateobs, throttle.(throt, observe.(sl))...)
        stateobs[] = [s[] for s in sl]
        sl
    end

    make_varchecks(model) = begin
        checks = PlotNested.plotchecks(model.records)
        map!((x...) -> [x...], varobs, observe.(checks)...)
        varobs[] = [c[] for c in checks]
        checks
    end

    make_envchecks(model) = begin
        checks = PlotNested.plotchecks(model.environment)
        map!((x...) -> [x...], envobs, observe.(checks)...)
        envobs[] = [c[] for c in checks]
        checks
    end

    flux_grids = (make_grid(STATE, TRANS), make_grid(STATE1, TRANS1,),
                  make_grid(STATE, TRANS), make_grid(STATE1, TRANS1));
    fluxbox = hbox(arrange_grid.(flux_grids)...)
    fluxobs = map((g...) -> g, observe_grid.(flux_grids)...)

    tstopobs = map(throttle(throt, observe(tstoptext))) do t
        int = tryparse(Int, t)
        int == nothing ? 8760hr : int * hr
    end
    envstartobs = throttle(throt, observe(envstart))

    map!(make_paramsliders, paramsliders, modelobs)
    map!(make_varchecks, varchecks, modelobs)
    map!(make_envchecks, envchecks, modelobs)
    map!(make_statesliders, statesliders, modelobs)

    map!(c -> vbox(hbox(c...)), varbox, varchecks)
    halfenv = min(12, length(envchecks[]))
    map!(c -> vbox(hbox(c[1:halfenv]...), hbox(c[halfenv+1:end]...)), envbox, envchecks)
    map!(s -> spreadwidgets(s), parambox, paramsliders)
    map!(s -> spreadwidgets(s), statebox, statesliders)

    map!(make_plot, plotobs, modelobs, stateobs[], solplotobs, varobs, envobs,
         fluxobs, observe(plottemp), observe(plotscale), observe(plotphoto), observe(plotpot), tstopobs[])
    map!(sol_plot, solplotobs, modelobs, paramobs, stateobs, tstopobs, envstartobs, app)

    modelobs[] = load_model(modeldrop[], envdrop[])

    ui = vbox(hbox(reload, dropbox), controlbox, fluxbox, varbox, envbox, plotobs, parambox, statebox);
    # dom"div"()
end

function colorforeach(parents)
    parentlist = union(parents)
    pallete = get(ColorSchemes.pastel, 0:(1/length(parentlist)):1)
    [pallete[indexin([parent], parentlist)[1]] for parent in parents]
end

state_slider(label, val) =
    slider(vcat(exp10.(range(-6, stop = 1, length = 100))) * mol, label=label, value=val)

make_grid(rownames, colnames) = begin
    rows = length(rownames)
    cols = length(colnames)
    [checkbox(false, label = join([string(rownames[r]), string(colnames[c])], ",")) for r = 1:rows, c = 1:cols]
end

arrange_grid(a) = hbox((vbox(a[:,col]...) for col in 1:size(a, 2))...)

observe_grid(a) = map((t...) -> [t[i + size(a,1) * (j - 1)] for i in 1:size(a,1), j in 1:size(a,2)], observe.(a)...)


electronapp(app; zoom=0.6) = begin
    ui = app(nothing)
    w = Window(Dict("webPreferences"=>Dict("zoomFactor"=>0.6)));
    # Blink.AtomShell.@dot w webContents.setZoomFactor($zoom)
    body!(w, ui);
end

webapp(app; port=8000) = webio_serve(page("/", req -> app(req)), port)


savecode(app, name) = begin
    lines = split("models[:$name] = " * codify(app.savedmodel), "\n")
    code = join([lines[1], "    environment = first(values(environments)),", lines[2:end]...], "\n")
    write("models/$name.jl", code)
end
