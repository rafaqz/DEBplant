dir = "DEBSCRIPTS" in keys(ENV) ? ENV["DEBSCRIPTS"] : pwd()

# Line Plots ####################################################

include(joinpath(dir, "src/plotmethods.jl"))

# Load environments at t1
environments, _ = loadenvironments(dir)
environment = environments[:t1]
envstart = round(typeof(1hr), STARTMONTH * MONTH_HOURS)

# Load models from source code
vars = (PlottableVars(), PlottableVars())
models = OrderedDict()
include.(readdir(joinpath(dir, "models"); join=true));
model = models[:bb];
u = init_state(model)
model = set_allometry(model, u);

# Set up plots
# theme(:sand)
# theme(:solarized)
# theme(:juno)
# theme(:solarized_light)
theme(:wong2)
gr()
# plotly()


# Single-simulation plots
plot_crossover(model, environment, u, envstart, 2)
savefig("plots/crossover")
plot_assim(model, environment, u, envstart, 6, 1:1)
savefig("plots/growth")
plot_assim(model, environment, u, envstart, 6, 1:2)
savefig("plots/growthassim")
plot_assim(model, environment, u, envstart, 6, 1:3)
savefig("plots/assimall")
plot_assim(model, environment, u, envstart, 1, 1:1)
savefig("plots/assimshort")
plot_assim(model, environment, u, envstart, 6, 3:5)
savefig("plots/assimswp")
plot_assim(model, environment, u, envstart, 6, 1:5)
savefig("plots/assimall")


# Multi-simulation plots
gr()
plot_years(model, environments, u, envstart)
savefig("plots/all")
plot_years(model, Dict(:t1=>environments[:t1]), u, 1.0hr)
savefig("plots/t1")
plot_years(model, Dict(:t2=>environments[:t2]), u, 1.0hr, "t2")
savefig("plots/t2")
plot_years(model, Dict(:t3=>environments[:t3]), u, 1.0hr, "t3")
savefig("plots/t3")


# Plot temperature response curve
temps = 0.0°C:0.1K:50.0°C
plot(x -> tempcorr(tempcorr_pars(model.shared), K(x)), temps;
    legend=false,
    ylabel="Correction",
    xlabel="Temperature"
)
savefig("plots/tempcorr")



# Map #################################################################

include(joinpath(dir, "src/mapping.jl"))

datapath = "/home/raf/Data/microclim"
radpath = joinpath(datapath, "SOLR/SOLR_2001.nc")
isdir(datapath) || error("Need to set datapath to you microclim dataset folder")
long = NCDatasets.Dataset(ds -> Array(ds["longitude"]), radpath)
lat = NCDatasets.Dataset(ds -> Array(ds["latitude"]), radpath)

skip = (:snowdepth, :soilwatercontent) 
MONTH_HOURS = 365.25 / 12 * 24hr
years = 2005:2010
shade = 0

# Import models
vars = (Vars(), Vars())
models = OrderedDict()
modeldir = joinpath(dir, "models")
include.(joinpath.(Ref(modeldir), readdir(modeldir)));
model = deepcopy(models[:bb]);
model.environment_start[] = oneunit(model.environment_start[])

# Load Australian border shapefile
shapefile = joinpath(dir, "data/ausborder_polyline.shp")
shp = open(shapefile) do io
    read(io, Shapefile.Handle)
end

if isfile(joinpath(dir, "data/yearly_outputs.jld2"))
    yearly_outputs = load("data/yearly_outputs.jld2", "yearly_outputs")
else
    @time yearly_outputs = map(y -> run_year(y, datapath, shade, model, skip), years)
    save("data/yearly_outputs.jld2", Dict("yearly_outputs" => yearly_outputs))
end


# Simple outline plots
points = (getindex.(Ref(long), [65, 60, 55]), getindex.(Ref(lat), [35, 35, 35]))
shape_plot(long, lat, points, ["T1", "T2", "T3"], (800,600), 3)
savefig("plots/shape.png")
nswlong = long[45:end]
nswlat = lat[30:46]
shape_plot(nswlong, nswlat, points, ["T1", "T2", "T3"], (300,300), 5)
savefig("plots/nsw.png")

# Maxium structural mass plot
gr()
maximum(skipmissing(yearly_outputs[1]))
year_sums = combine_year.(yearly_outputs)
plts = build_plot.(year_sums, string.(years), (false, true, false, true, false, true))
maps = plot(plts...; layout=grid(length(plts)÷2, 2, widths=[0.423, 0.577]), size=(1000,1300), dpi=100)
