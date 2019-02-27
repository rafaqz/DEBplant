using Revise, Unitful, Microclimate, JLD2, DataStructures, Flatten, FieldMetadata, OrdinaryDiffEq
using Photosynthesis, DynamicEnergyBudgets
using DataStructures

import Plots:px, pct, GridLayout
using Photosynthesis: potential_dependence
using DynamicEnergyBudgets: STATE, STATE1, TRANS, TRANS1, shape_correction, define_organs,
      photosynthesis, split_state, HasCN, HasCNE, has_reserves, tempcorr_pars,
      assimilation_pars, assimilation_vars, parconv, w_V, build_vars
using Unitful: °C, K, Pa, kPa, MPa, J, kJ, W, L, g, kg, cm, m, s, hr, d, mol, mmol, μmol, σ

using Microclimate, Unitful, DataStructures, Setfield, DataStructures, JLD2

const STATEKEYS = (:PS, :VS, :MS, :CS, :NS, :ES, :PR, :VR, :MR, :CR, :NR, :ER)
const STATELABELS = tuple(vcat([string("Shoot ", s) for s in STATE], [string("Root ", s) for s in STATE])...)

loadenvironments(dir) = begin
    locationspath = joinpath(dir, "microclimate/locations.jld")
    @load locationspath t1 t2 t3 t4
    environments = OrderedDict(:t1 => t1, :t2 => t2, :t3 => t3)
    env = t1
    tspan = (0:1:length(radiation(env)) - 1) * hr
    environments, tspan
end

# The zero crossing of allometry is the seed size.
set_allometry(model, state) = begin
    if :β0 in fieldnames(typeof(model.params[1].allometry_pars))
        model = @set model.params[1].allometry_pars.β0 = state[2] * w_V(model.shared) * 0.999999
    end
    if :β0 in fieldnames(typeof(model.params[2].allometry_pars))
        model = @set model.params[2].allometry_pars.β0 = state[8] * w_V(model.shared) * 0.999999
    end
    model
end

assimvars(::AbstractCAssim) = DynamicEnergyBudgets.ShootVars()
assimvars(::AbstractNAssim) = DynamicEnergyBudgets.RootVars()
assimvars(::FvCBPhotosynthesis) = DynamicEnergyBudgets.FvCBShootVars()

function update_vars(m, tstop)
    m = @set m.records[1].vars = build_vars(assimvars(m.params[1].assimilation_pars), 1:ustrip(tstop))
    m = @set m.records[2].vars = build_vars(assimvars(m.params[2].assimilation_pars), 1:ustrip(tstop))
    m
end
