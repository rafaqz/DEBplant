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

# Set up the plot
plotly() # Use Plotly for plotting - can be zoomed in
# gr() # Use GR for plotting
# theme(:sand)
# theme(:solarized)
# theme(:juno)
theme(:solarized_light)
# theme(:wong2)
plotsize = (2000, 1000) # Set plot size in ui

# Build the app interface
app = ModelApp(models, environments, plotsize, tspan, nothing);

# Run the app in Atom
# display(app)

# Or create an electron desktop app
electronapp(app; zoom=0.4)
