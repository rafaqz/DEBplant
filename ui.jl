# Use a set script pat otherwise the current working directory
dir = "DEBSCRIPTS" in keys(ENV) ? ENV["DEBSCRIPTS"] : pwd()
# Load the app scripts
include(joinpath(dir, "app.jl"))
 # Import environments 
environments, tspan = loadenvironments(dir);
environments[:controls] = MicroclimControl();
# vars = (Vars(), Vars())
vars = (PlottableVars(), PlottableVars())

FieldMetadata.flattenable(::AbstractMicroclimate, x) = false
FieldMetadata.flattenable(::AbstractMicroclimControl, x) = true

# u = [1e-4, 1e-4, 1.0, 1e-4, 1e-4, 10.0]u"mol" # Initial value
du = fill(0.0u"mol/hr", 6)

# Import models
models = OrderedDict()
modeldir = joinpath(dir, "models")
include.(joinpath.(Ref(modeldir), readdir(modeldir)));

plant = models[:bb];
@set! plant.environment = environments[:t1];
plant.records[1].vars

# plant(du, u, nothing, 1u"hr")
# plant.dead[] = false
# prob = DiscreteProblem(plant, ustrip(u), (0, 8000));
# sol = solve(prob, FunctionMap(scale_by_time = true));
# plot(sol)

# Build the app interface
plotsize = (2000, 1000)
plotly()
app = ModelApp(models, environments, plotsize, tspan, nothing);

# Electron desktop app
electronapp(app; zoom=0.6)


# Save the code for the current state of the model in the app 
# savecode(app, "maturity")
