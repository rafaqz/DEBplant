using Revise, Setfield
using Unitful, Microclimate, DataStructures, Flatten, FieldMetadata, 
      LabelledArrays, OrdinaryDiffEq, Photosynthesis, DynamicEnergyBudgets, 
      DataStructures, Plots, ColorSchemes, UnitfulPlots, JLD2

import Plots:px, pct, GridLayout
using DynamicEnergyBudgets: STATE, STATE1, TRANS, TRANS1, shape_correction, define_organs,
      photosynthesis, split_state, HasCN, HasCNE, has_reserves, tempcorr_pars,
      assimilation_pars, parconv, w_V, build_vars, allometry_pars
using Unitful: °C, K, Pa, kPa, MPa, J, kJ, W, L, g, kg, g, mg, cm, m, s, hr, d, mol, mmol, μmol, σ

# const STATEKEYS = (:PS, :VS, :MS, :CS, :NS, :ES, :PR, :VR, :MR, :CR, :NR, :ER)
const STATEKEYS = tuple([Symbol(string(s, :S)) for s in STATE]..., [Symbol(string(s, :R)) for s in STATE]...)
const STATELABELS = tuple([string("Shoot ", s) for s in STATE]..., [string("Root ", s) for s in STATE]...)

loadenvironments(dir) = begin
    locationspath = joinpath(dir, "microclimate/locations.jld")
    @load locationspath t1 t2 t3
    environments = OrderedDict{Symbol,Any}(:t1 => t1, :t2 => t2, :t3 => t3)
    env = t1
    tspan = (0:1:length(radiation(env)) - 1) * hr
    environments, tspan
end

# The zero crossing of allometry is the seed size.
set_allometry(model, state) = begin
    if :β0 in fieldnames(typeof(model.params[1].allometry_pars))
        model = @set model.params[1].allometry_pars.β0 = state[:VS] * w_V(model.shared) * 0.999999
    end
    if :β0 in fieldnames(typeof(model.params[2].allometry_pars))
        model = @set model.params[2].allometry_pars.β0 = state[:VR] * w_V(model.shared) * 0.999999
    end
    model
end


function update_vars(m, tstop)
    m = @set m.records[1].vars = build_vars(Vars(), 1:ustrip(tstop))
    m = @set m.records[2].vars = build_vars(Vars(), 1:ustrip(tstop))
    m
end

function discrete_solve(model, u0, tstop)
    u = deepcopy(u0)
    mx = deepcopy(u0)
    du = u ./ unit(tstop)
    for i = 2oneunit(tstop):oneunit(tstop):tstop
        model(du, u, nothing, i)
        model.dead[] && break
        u .+= du .* unit(tstop)
        mx = max.(mx, u)
    end
    # Return the maximum values
    mx
end

# function discrete_solve(model, u, tstop)
    # prob = DiscreteProblem(model, u, (oneunit(tstop), tstop))
    # solve(prob, FunctionMap(scale_by_time = true))
# end
