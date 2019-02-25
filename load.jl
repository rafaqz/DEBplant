using Microclimate, Unitful, DataStructures, Setfield, DataStructures, JLD2

const STATEKEYS = (:PS, :VS, :MS, :CS, :NS, :ES, :PR, :VR, :MR, :CR, :NR, :ER)
const STATELABELS = tuple(vcat([string("Shoot ", s) for s in STATE], [string("Root ", s) for s in STATE])...)

loadenvironments(dir) = begin
    locationspath = joinpath(dir, "microclimate/locations.jld")
    @load locationspath t1 t2 t3 t4
    environments = OrderedDict(:t1 => t1, :t2 => t2, :t3 => t3, :t4 => t4)
    env = t1
    tspan = (0:1:length(radiation(env)) - 1) * hr
    environments, tspan
end

# The zero crossing of allometry is the seed size.
set_allometry(model, state) = begin
    @set! model.params[1].allometry_pars.β0 = state[2] * w_V(model.shared)
    @set! model.params[2].allometry_pars.β0 = state[8] * w_V(model.shared)
    model
end
