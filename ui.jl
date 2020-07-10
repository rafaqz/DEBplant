#=
A DEB model user interface. This allows manual construction
and parametrisation of a DEB plant model, run with forcing variables
from microclimate data.

You can swap out any model components by selecting them from the dropdowns
and hitting the "Reload" button. It should give you an updated set of sliders
to control the new model with.

Giving the model a name and hitting save will save the load script for the model in the
models folder. This will be available next time you load the ui.
If you use an existing name, it will be written over.

If the plot is white initially, just click a flux or vars 
check box to make it build.

WARNING: this will take a while to load.
=#

dir = "DEBSCRIPTS" in keys(ENV) ? ENV["DEBSCRIPTS"] : pwd()

# Load the app scripts
include(joinpath(dir, "src/app.jl"))

 # Import environments 
environments, tspan = loadenvironments(dir);
environments[:controls] = MicroclimControl();

# Import models
vars = (PlottableVars(), PlottableVars())
models = OrderedDict()
modeldir = joinpath(dir, "models")
include.(joinpath.(Ref(modeldir), readdir(modeldir)));

# Build the app interface
plotly() # Use Plotly.js. Alternatively use gr()
plotsize = (2000, 1000) # Set plot size in ui
app = ModelApp(models, environments, plotsize, tspan, nothing);

# Create an electron desktop app
electronapp(app; zoom=0.6)
