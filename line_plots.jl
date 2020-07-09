dir = "DEBSCRIPTS" in keys(ENV) ? ENV["DEBSCRIPTS"] : pwd()

include("multiplot.jl")

# Load environments at t1
environments, _ = loadenvironments(dir)
environment = environments[:t1]
envstart = round(typeof(1hr), STARTMONTH * MONTH_HOURS)

# Load models
vars = (PlottableVars(), PlottableVars())
models = OrderedDict()
include.(readdir(joinpath(dir, "models"); join=true));
model = models[:bb];
u = init_state(model)
model = set_allometry(model, u);

# Choose a plot theme
# theme(:sand)
# theme(:solarized)
# theme(:juno)
# theme(:solarized_light)
theme(:wong2)

# savefig("plots/tempcorr")

# Single-simulation plots
plot_crossover(model, environment, u, envstart, 2)
plot_assim(model, environment, u, envstart, 6, 1:1, "growth")
plot_assim(model, environment, u, envstart, 6, 1:2, "growthassim")
plot_assim(model, environment, u, envstart, 6, 1:3, "assimall")
plot_assim(model, environment, u, envstart, 1, 1:1, "assimshort")
plot_assim(model, environment, u, envstart, 6, 3:5, "assimswp")
plot_assim(model, environment, u, envstart, 6, 1:5, "assimall")


# Multi-simulation plots
gr()
model = set_allometry(models[:bb], u)
plot_years(model, environments, u, envstart, "all")
# plot_years(model, Dict(:t1=>environments[:t1]), u, 1.0hr, "t1")
# plot_years(model, Dict(:t2=>environments[:t2]), u, 1.0hr, "t2")
# plot_years(model, Dict(:t3=>environments[:t3]), u, 1.0hr, "t3")


# Plot temperature response curve
temps = 0.0°C:0.1K:50.0°C
plot(x -> tempcorr(tempcorr_pars(model.shared), K(x)), temps;
    legend=false,
    ylabel="Correction",
    xlabel="Temperature"
)
