using Revise
using Unitful
using UnitlessFlatten
using OrdinaryDiffEq
using DynamicEnergyBudgets
using Photosynthesis
using AxisArrays
using InteractBulma, InteractBase, Blink, WebIO, Observables, CSSUtil
using Mux
using FieldMetadata
# using IterableTables, DataFrames, TypedTables
# using JLD2
using Microclimate
using Plots, UnitfulPlots, PlotNested, StatPlots
import Plots:px, pct, GridLayout
import InteractBase: WidgetTheme, libraries
using DynamicEnergyBudgets: STATE, STATE1, TRANS, TRANS1, scaling, define_organs, photosynthesis, split_state

dir = "/home/raf/julia/DynamicEnergyBudgets/scratch/"
# dir = "/home/cloud"

# struct MyTheme <: WidgetTheme; end
# libraries(::MyTheme) = vcat(libraries(Bulma()), [joinpath(dir, "custom.css")])
# settheme!(MyTheme())
# gr()
pyplot()

# function load_environment()
#     # using IndexedTables
#     environment = load(joinpath(dir, "environment.jld"))["environment"]
#     # environment = nichemap_global("Adelaide", years=10)
#     env2 = Microclimate.MicroclimateTable(
#       Table(environment.soil),
#       Table(environment.shadsoil),
#       Table(environment.metout),
#       Table(environment.shadmet),
#       Table(environment.soilmoist),
#       Table(environment.shadmoist),
#       Table(environment.humid),
#       Table(environment.shadhumid),
#       Table(environment.soilpot),
#       Table(environment.shadpot),
#       Table(environment.plant),
#       Table(environment.shadplant),
#       environment.RAINFALL,
#       environment.dim,
#       environment.ALTT,
#       environment.REFL,
#       environment.MAXSHADES,
#       environment.longlat,
#       environment.nyears,
#       environment.timeinterval,
#       environment.minshade,
#       environment.maxshade,
#       environment.DEP,
#     )
# end

const t = 0:1:10000
const statelabels = vcat([string("Shoot ", s) for s in STATE], [string("Root ", s) for s in STATE])
const u12 = [0.0, 1e-4, 0.0, 1e-4, 1e-4, 1e-4, 0.0, 1e-4, 0.0, 1e-4, 1e-4, 10.0]
const u18 = [0.0, 1e-4, 0.0, 1e-4, 1e-4, 1e-4, 0.0, 1e-4, 0.0, 1e-4, 1e-4, 1e-4, 0.0, 1e-4, 0.0, 1e-4, 1e-4, 10.0]

photo(o, u, x) = begin
    update_photovars(o.vars.assimilation, x)
    ustrip(photosynthesis(o.params.assimilation, o, u))
end

update_photovars(v::CarbonVars, x) = v.J_L_F = x * oneunit(v.J_L_F)
update_photovars(v::Photosynthesis.PhotoVars, x) = v.par = x * oneunit(v.par)

function make_plot(model, checkobs, sliderobs, fluxobs, timespan)
    # namedparams = AxisArray([params...], Axis{:parameters}(names))
    # try
    length(sliderobs) > 0 || return plot()
        m = reconstruct(model, sliderobs)
        organs = define_organs(m, 1)
        o = organs[1]
        state = split_state(organs, u12 .* 10000, 0)
        dataplots = plot_selected(model, checkobs, 1:timespan)
        prob = DiscreteProblem(m, u12, (1, timespan))
        sol = solve(prob, FunctionMap(scale_by_time = true))
        # solplot1 = plot(sol, vars = [1:6...], plotdensity=400, legend=:topleft, labels=statelabels[1:6], ylabel="State (CMol)", 
                      # xlabel=string(model.params[1].name, " : ",  typeof(model.params[1].assimilation).name, " - time (hr)"))
        # solplot2 = plot(sol, vars = [7:12...], plotdensity=400, legend=:topleft, labels=statelabels[7:12], ylabel="State (CMol)", 
                      # xlabel=string(model.params[2].name, " : ", typeof(model.params[2].assimilation).name, " - time (hr)"))
        solplot1 = plot(sol.t, sol', legend=:topleft, labels=reshape(statelabels, 1, length(statelabels)), ylabel="State (CMol)", 
                      xlabel=string(model.params[1].name, " : ",  typeof(model.params[1].assimilation).name, " - time (hr)"))
        # solplot2 = plot(sol.t, sol', vars = [7:12...], plotdensity=400, legend=:topleft, labels=statelabels[7:12], ylabel="State (CMol)", 
                      # xlabel=string(model.params[2].name, " : ", typeof(model.params[2].assimilation).name, " - time (hr)"))
        solplot = plot(solplot1, dataplots..., layout=Plots.GridLayout(1+length(dataplots), 1))
        tempcorr(10u"°C", m.shared.tempcorr)
        arrh = plot(x -> ustrip(tempcorr(x * 1.0u"°C", m.shared.tempcorr)), 0.0:1.0:50.0, legend=false, ylabel="Correction", xlabel="°C")
        scaleplot(i) = plot(x -> ustrip(scaling(m.params[i].scaling, x * 1.0u"mol")), 0.0:0.01:100.0, legend=false, ylabel="Correction", xlabel="CMols Stucture")
        scaling1 = scaleplot(1)
        scaling2 = scaleplot(2)
        photoplot = plot(x -> photo(organs[1], state[1], x), 0:10:4000) 
        funcs = plot(arrh, scaling1, scaling2, photoplot, layout=Plots.GridLayout(4, 1))
        l = Plots.GridLayout(1, 2)
        l[1, 1] = GridLayout(1, 1, width=0.8pct)
        l[1, 2] = GridLayout(1, 1, width=0.2pct)
        plot(solplot, funcs, size=(1700, 800), layout=l) #, fmt = :png
    # catch err
        # solplot = plot(xlabel="model run failed", size=(1700, 800))
    # end
end

function spreadwidgets(widgets; cols = 7)
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

function muxapp(req) # an "App" takes a request, returns the output
    env2 = nothing #load_environment()

    noenv = button("No Environment")
    plantcne = button("PlantCNE")
    plantcn = button("PlantCN")
    no_maturity = button("No Maturity")
    constant = button("Constant Assimilation")
    fvcbplant = button("Farquhar")
    # plant3 = button("3 organ plant")
    # fvcbplant3 = button("3 organ Farquhar")
    timespan = slider(20:9999, label="Timespan")

    model = Observable{Any}(DynamicEnergyBudgets.Plant(environment=env2, time=t))
    u = Observable{Vector{Float64}}(u12)
    plt = Observable{Any}(0);
    sliders = Observable{Vector{Widget{:slider}}}(Widget{:slider}[]);
    sliderobs = Observable{Vector{Float64}}(Float64[]);
    sliderbox = Observable{Any}(dom"div"());
    checks = Observable{Vector{Widget{:checkbox}}}(Widget{:checkbox}[]);
    checkobs = Observable{Vector{Bool} where N}(Bool[]);
    checkboxbox = Observable{Any}(0);


    # map!(x -> u12, u, observe(plant));
    # map!(x -> u12, u, observe(constant));
    # map!(x -> u12, u, observe(no_maturity));
    # map!(x -> u12, u, observe(fvcbplant));
    # map!(x -> u18, u, observe(fvcbplant3));
    # map!(x -> u18, u, observe(plant3));

    map!(x -> DynamicEnergyBudgets.Plant(time=t), model, observe(noenv));
    map!(x -> DynamicEnergyBudgets.Plant(environment=env2, time=t), model, observe(plantcne));
    map!(x -> DynamicEnergyBudgets.PlantCN(environment=env2, time=t), model, observe(plantcn));
    map!(x -> DynamicEnergyBudgets.ConstantPlant(environment=env2, time=t), model, observe(constant));
    map!(x -> DynamicEnergyBudgets.NoMaturityPlant(environment=env2, time=t), model, observe(no_maturity));
    map!(x -> DynamicEnergyBudgets.FvCBPlant(environment=env2, time=t), model, observe(fvcbplant));
    # map!(x -> DynamicEnergyBudgets.Plant3(environment=env2, time=t), model, observe(plant3));
    # map!(x -> DynamicEnergyBudgets.FvCBPlant3(environment=env2, time=t), model, observe(fvcbplant3));

    function make_checks(m)
        checks = plotchecks(m)
        map!((x...) -> [x...], checkobs, observe.(checks)...)
        checks
    end
    map!(make_checks, checks, model);

    flux_grids = (make_grid(STATE, TRANS), make_grid(STATE1, TRANS1,),
                  make_grid(STATE, TRANS), make_grid(STATE1, TRANS1));
    flux_interface = hbox(arrange_grid.(flux_grids)...)
    fluxobs = map((g...) -> g, observe_grid.(flux_grids)...)

    function make_sliders(m)
        params = flatten(Vector, m)
        fnames = fieldnameflatten(Vector, m)
        parents = parentflatten(Vector, m)
        limits = metaflatten(Vector, m, DynamicEnergyBudgets.limits)
        descriptions = metaflatten(Vector, m, DynamicEnergyBudgets.description)
        attributes = broadcast((p, n, d) -> Dict(:title => "$p.$n: $d"), parents, fnames, descriptions)
        sl = broadcast((x,l,v,a) -> InteractBase.slider(x[1]:(x[2]-x[1])/100:x[2], label=string(l), value=v, attributes=a), 
                       limits, fnames, params, attributes)
        # description = metaflatten(Vector, m, DynamicEnergyBudgets.label)
        # InteractBase.tooltip!.(sliders, description) 
        map!((x...) -> [x...], sliderobs, throttle.(0.1, observe.(sl))...) 
        sl
    end

    map!(make_sliders, sliders, model);
    map!(make_plot, plt, model, checkobs, sliderobs, fluxobs, throttle(0.1, observe(timespan)))

    map!(c -> hbox(c...), checkboxbox, checks)
    map!(s -> spreadwidgets(s), sliderbox, sliders)
    model[] = DynamicEnergyBudgets.Plant(time=t);
    o = observe(sliders[][5])
    o[] = 0.45699 
    modelbox = hbox("Model: ", noenv, plantcne, plantcn, no_maturity, constant, fvcbplant, timespan)
    ui = vbox(plt, modelbox, sliderbox, checkboxbox, flux_interface);
end

make_grid(rownames, colnames) = begin
    rows = length(rownames)
    cols = length(colnames)
    [checkbox(false, label = join([string(rownames[r]), string(colnames[c])], ",")) for r = 1:rows, c = 1:cols]
end

arrange_grid(a) = hbox((vbox(a[:,col]...) for col in 1:size(a, 2))...)

observe_grid(a) = map((t...) -> [t[i + size(a,1) * (j - 1)] for i in 1:size(a,1), j in 1:size(a,2)], observe.(a)...)

# Compile everything first
env2 = nothing # load_environment()
a = DynamicEnergyBudgets.ConstantPlant(time=t)
a = DynamicEnergyBudgets.ConstantPlant(environment=env2, time=t)
b = DynamicEnergyBudgets.NoMaturityPlant(environment=env2, time=t)
c = DynamicEnergyBudgets.FvCBPlant(environment=env2, time=t)
d = DynamicEnergyBudgets.Plant(environment=env2, time=t)
# e = DynamicEnergyBudgets.Plant3(environment=env2, time=t)
# f = DynamicEnergyBudgets.FvCBPlant3(environment=env2, time=t)
make_plot(a, (true), flatten(a), (), 1000)
make_plot(b, (true), flatten(b), (), 1000)
make_plot(c, (true), flatten(c), (), 1000)
make_plot(d, (true), flatten(d), (), 1000)
# make_plot(e, u18, (true), flatten(e), 1000)
# make_plot(f, u18, (true), flatten(f), 1000)
# muxapp(1)

w = Window(Dict("webPreferences"=>Dict("zoomFactor"=>0.7)));
Blink.AtomShell.@dot w webContents.setZoomFactor(0.7)
ui = muxapp(nothing) 
body!(w, ui);
# opentools(w)

# webio_serve(page("/", req -> muxapp(req)), 8000)

# function accordian(sliderbox)
#  dom"div"(
#   dom"section[class=accordions]"(
#       dom"article[class=accordion is-active]"(
#           dom"div[class=accordion-header toggle]"( "Sliders"),
#           dom"div[class=accordion-body]"(
#              dom"div[class=accordion-content]"(sliderbox))),
#       dom"article[class=accordion]"(
#           dom"div[class=accordion-header]"("Plots"),
#           dom"div[class=accordion-body]"(
#               dom"div[class=accordion-content]"(
#                  "content")))),
#     )
# end

# @ji w "var accordions = bulmaAccordion.attach();"
