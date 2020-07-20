
#= Map generation #################################################################

Saved .jld2 files are included for generating maps without running simulations.
This will happen by default.

To run simulations you will need to get the MicroclimOz dataset,
here using the commented out `download_microclim` command. Then
set `run_sims = true`.

Downloading files will take quite a few hours, and requires 380GB of storage 
when unzipped. After that, simulations may take an hour or more to run.
=#

dir = "DEBSCRIPTS" in keys(ENV) ? ENV["DEBSCRIPTS"] : pwd()
include(joinpath(dir, "src/mapping.jl"))


# Load environment -----------------------------------------------------

# Download data. 
# datapath = "/home/raf/Data/microclim_oz"
# download_microclim(datapath; overwrite=false)

# Import spatial data

# From netcdf files
# radpath = joinpath(datapath, "SOLR/SOLR_2001.nc")
# lons = NCDatasets.Dataset(ds -> Array(ds["longitude"]), radpath)
# lats = reverse(NCDatasets.Dataset(ds -> Array(ds["latitude"]), radpath))
# save(joinpath(dir, "data/coords.jld2"), Dict("lats" => lats, "lons" => lons))

# From saved jld2
lats, lons = load(joinpath(dir, "data/coords.jld2"), "lats", "lons")


# Load Australian border shapefile
shapefile = joinpath(dir, "data/ausborder_polyline.shp")
shp = open(io -> read(io, Shapefile.Handle), shapefile)

# Transect locations
transect_lons = getindex.(Ref(lons), [65, 60, 55])
transect_lats = getindex.(Ref(lats), [20, 20, 20])


# Import models --------------------------------------------------------

environment = nothing
vars = (Vars(), Vars())
models = OrderedDict()
modeldir = joinpath(dir, "models")
include.(joinpath.(Ref(modeldir), readdir(modeldir)));
model = deepcopy(models[:bb]);
model.environment_start[] = oneunit(model.environment_start[])


# Simulations ----------------------------------------------------------

run_sims = false
skip = (:snowdepth, :soilwatercontent) 
MONTH_HOURS = 365.25 / 12 * 24hr
years = 2005:2010
shade = 0

# Run simulations or load previous run from disk 
if !run_sims && isfile(joinpath(dir, "data/yearly_outputs.jld2"))
    yearly_outputs = load(joinpath(dir, "data/yearly_outputs.jld2"), "yearly_outputs")
else
    isdir(datapath) || error("Need to set `datapath` to you microclim dataset folder")
    @time yearly_outputs = map(y -> run_year(y, datapath, shade, model, skip), years)
    save(joinpath(dir, "data/yearly_outputs.jld2"), Dict("yearly_outputs" => yearly_outputs))
end

# Simple map outlines -------------------------------------------------

shape_plot(lons, lats, transect_lons, transect_lats, ["T1", "T2", "T3"], (800,600), 3)
savefig("plots/shape.png")
nswlon = lons[45:end]
nswlat = lats[9:24]
shape_plot(nswlon, nswlat, transect_lons, transect_lats, ["T1", "T2", "T3"], (300,300), 5)
savefig(joinpath(dir, "plots/nsw.png"))


# Maximum structural mass map -----------------------------------------

dir = "DEBSCRIPTS" in keys(ENV) ? ENV["DEBSCRIPTS"] : pwd()
include(joinpath(dir, "src/mapping.jl"))
month_sums = combine_month(yearly_outputs)
year_sums = combine_year.(yearly_outputs)

# Monthly breakdown
month_legends = (false, false, false, false, false, false, false, false, false, false, false, false) 
month_plts = map(month_sums, monthname.(1:12), month_legends) do sum, year, haslabel
    build_map(sum, year, haslabel, lons, lats, transect_lons, transect_lats)
end
month_map = plot(month_plts...; 
    layout=grid(4, 3), 
    size=(1000,1300), 
    dpi=100,
)
savefig(joinpath(dir, "plots/monthly_map.png"))

# Yearly breakdown
year_plts = map(year_sums, string.(years), (false, true, false, true, false, true)) do sum, year, haslabel
    build_map(sum, year, haslabel, lons, lats, transect_lons, transect_lats)
end
year_map = plot(year_plts...; 
    layout=grid(length(year_plts) ÷ 2, 2; widths=[0.423, 0.577]), 
    size=(1000,1300), 
    dpi=100
)
savefig(joinpath(dir, "plots/yearly_map.png"))
