dir = "DEBSCRIPTS" in keys(ENV) ? ENV["DEBSCRIPTS"] : pwd()
include(joinpath(dir, "app.jl"))

statelabels = tuple(vcat([string("Shoot ", s) for s in STATE], [string("Root ", s) for s in STATE])...)

# Import environments 
locationspath = joinpath(dir, "microclimate/locations.jld")
@load locationspath tas desert qld
environments = OrderedDict(:Tas=>tas, :Desert=>desert, :QLD=>qld)
env = first(values(environments))

tspan = (0:1:length(radiation(env)) - 1) * hr

# Load all the saved models
models = OrderedDict()
modeldir = joinpath(dir, "models")
include.(joinpath.(Ref(modeldir), readdir(modeldir)));
models

app = ModelApp(models, environments, tspan, statelabels, nothing);

# Electron desktop app
electronapp(app; zoom=0.6)

# Save the code for the current state of the model in the app 
# It will be loaded automatically the next time this script is run
# Choose a new name if this is a different model, leave if it's an update
# savecode(app, "maturity")
