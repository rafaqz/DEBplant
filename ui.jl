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

# Build the app interface
app = ModelApp(models, environments, tspan, nothing);

# Electron desktop app
electronapp(app; zoom=0.5)



# Save the code for the current state of the model in the app 
# It will be loaded automatically the next time this script is run
# Choose a new name if this is a different model, leave if it's an update
# savecode(app, "maturity")
