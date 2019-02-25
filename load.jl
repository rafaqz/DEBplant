using Microclimate, Unitful, DataStructures, Setfield, DataStructures, JLD2

loadenvironments(dir) = begin
    locationspath = joinpath(dir, "microclimate/locations.jld")
    @load locationspath environments
    environments
    env = first(values(environments))
    tspan = (0:1:length(radiation(env)) - 1) * hr
    environments, tspan
end

# The zero crossing of allometry is the seed size.
set_allometry(model, state) = begin
    @set! model.params[1].allometry_pars.β0 = state.VS * w_V(model.shared)
    @set! model.params[2].allometry_pars.β0 = state.VR * w_V(model.shared)
end
