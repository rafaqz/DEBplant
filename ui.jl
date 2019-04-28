# Use a set script path otherwise the current working directory
dir = "DEBSCRIPTS" in keys(ENV) ? ENV["DEBSCRIPTS"] : pwd()
# Load the app scripts
include(joinpath(dir, "app.jl"))
 # Import environments 
environments, tspan = loadenvironments(dir)
environments[:controls] = MicroclimControl()

# Import models
models = OrderedDict()
modeldir = joinpath(dir, "models")
include.(joinpath.(Ref(modeldir), readdir(modeldir)));
plotsize = (2000, 1000)

# Build the app interface
app = ModelApp(models, environments, plotsize, tspan, nothing);

# Electron desktop app
electronapp(app; zoom=0.5)


# Save the code for the current state of the model in the app 
# savecode(app, "maturity")
