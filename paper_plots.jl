# Line Plots ####################################################

# Load helper scripts
dir = "DEBSCRIPTS" in keys(ENV) ? ENV["DEBSCRIPTS"] : pwd()
include(joinpath(dir, "src/plotmethods.jl"))

# Load environments at t1 ------------------------------------

# transect, _ = transect_from_netcdf(dir)
transect, _ = transect_from_saved(dir)
environment = transect[:t1]
envstart = round(typeof(1hr), STARTMONTH * MONTH_HOURS)

# Load models from source code ------------------------------

#= These are added to the models Dict when scripts are included. 
It's a little clunky, but it lets us use models saved from the UI. =#
vars = (PlottableVars(), PlottableVars())
models = OrderedDict()
include.(readdir(joinpath(dir, "models"); join=true));
model = models[:bb];
u = init_state(model)
model = set_allometry(model, u);


# Plots -------------------------------------------------------
gr()

# theme(:sand)
# theme(:solarized)
# theme(:juno)
# theme(:solarized_light)
theme(:default)
theme(:wong2)
theme(:vibrant)
theme(:mute)
theme(:bright)
theme(:wong2)

# Plot and save single-simulation plots
plot_crossover(model, environment, u, envstart, 2)
savefig("plots/crossover")
plot_byname(model, environment, u, envstart, 6, [:mass, :assim])
savefig("plots/growthassim")

# Plot and save multi-simulation plots
plot_years(model, transect, u, 1.0hr)
savefig("plots/transect_multiplot")
plot_years(model, Dict(:t1=>transect[:t1]), u, 1.0hr)
savefig("plots/t1_multiplot")
plot_years(model, Dict(:t2=>transect[:t2]), u, 1.0hr)
savefig("plots/t2_multiplot")
plot_years(model, Dict(:t3=>transect[:t3]), u, 1.0hr)
savefig("plots/t3_multiplot")

# Plot and save temperature response curve
temps = 0.0°C:0.1K:50.0°C
plot(x -> tempcorr(tempcorr_pars(model.shared), K(x)), temps;
    legend=false,
    ylabel="Correction",
    xlabel="Temperature"
)
savefig("plots/tempcorr")

