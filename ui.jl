# Use a set script pat otherwise the current working directory

dir = "DEBSCRIPTS" in keys(ENV) ? ENV["DEBSCRIPTS"] : pwd()

# Load the app scripts
include(joinpath(dir, "app.jl"))

 # Import environments 
environments, tspan = loadenvironments(dir);
environments[:controls] = MicroclimControl();

# Import models
vars = (PlottableVars(), PlottableVars())
models = OrderedDict()
modeldir = joinpath(dir, "models")
include.(joinpath.(Ref(modeldir), readdir(modeldir)));

plant = models[:bb];
@set! plant.environment = environments[:t1];

# Build the app interface
plotly() # Use Plotly.js. Alternatively use gr()
plotsize = (2000, 1000) # Set plot size in ui
app = ModelApp(models, environments, plotsize, tspan, nothing);

# Create an electron desktop app
electronapp(app; zoom=0.6)

# Save the code for the current state of the model in the app 
# savecode(app, "maturity")
