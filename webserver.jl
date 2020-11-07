# Settings
dir = dirname(@__FILE__)
port=8083
plotsize = (2000, 1000) # Set plot size in ui
# ui_theme = :sand
# ui_theme = :solarized
# ui_theme = :juno
# ui_theme = :solarized_light
# ui_theme = :wong2
ui_theme = :default


# Build the ui
println("Building ui...")
include(joinpath(dir, "src/app.jl"))

environments, tspan = transect_from_saved(dir);
environments[:controls] = MicroclimControl();
environment = first(environments)
environments

# Load models
vars = (PlottableVars(), PlottableVars())
models = OrderedDict()
include.(readdir(joinpath(dir, "models"); join=true));

# Set up the plot
plotly() # Use Plotly for plotting - can be zoomed in
# gr() # Use GR for plotting
theme(ui_theme)

# Start a server on port 80
println("Startinge server on port $port")

app = ModelApp(models, environments, plotsize, tspan, nothing);
webapp(app; port=port)
