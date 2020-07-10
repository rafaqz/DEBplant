using Setfield
using Unitful, Plots, UnitfulRecipes, Microclimate, DataStructures, Flatten, FieldMetadata, 
      LabelledArrays, OrdinaryDiffEq, Photosynthesis, DynamicEnergyBudgets, 
      DataStructures, ColorSchemes, JLD2, DimensionalData

using DynamicEnergyBudgets: scaling_correction, define_organs,
      photosynthesis, split_state, HasCN, HasCNE, has_reserves, tempcorr_pars,
      assimilation_pars, parconv, w_V, build_vars, allometry_pars

using Unitful: °C, K, Pa, kPa, MPa, J, kJ, W, L, g, kg, g, mg, cm, m, s, hr, d, mol, mmol, μmol, σ

import Plots:px, pct, GridLayout

loadenvironments(dir) = begin
    locationspath = joinpath(dir, "data/locations.jld")
    @load locationspath t1 t2 t3
    environments = OrderedDict{Symbol,Any}(:t1 => t1, :t2 => t2, :t3 => t3)
    env = t1
    tspan = (0:1:length(radiation(env)) - 1) * u"hr"
    environments, tspan
end

init_state(model::AbstractOrganism) = init_state(has_reserves.(define_organs(model, 1hr)), model)
init_state(::NTuple{2,HasCN}, model) = begin
    xdim = dims(first(model.records).J, X)
    xval = DimensionalData.unwrap(val(xdim))
    A = zeros(length(xdim) * length(model.records)) * mol
    newxval = Val((map(x -> Symbol(x, :S), xval)..., map(x -> Symbol(x, :R), xval)...))
    newxdim = X(newxval, mode(xdim), nothing)
    u = DimensionalArray(A, (newxdim,))
    u[:VS] = 0.2mg / (25.0g/mol)
    u[:CS] = 5.0mg  / (25.0g/mol)
    u[:NS] = 0.2mg  / (25.0g/mol)
    u[:VR] = 0.04mg / (25.0g/mol)
    u[:CR] = 1.0mg  / (25.0g/mol)
    u[:NR] = 0.04mg  / (25.0g/mol)
    u
end
init_state(::NTuple{2,HasCNE}, model) = [0.0, 1e-4, 0.0, 1e-4, 1e-4, 1e-4, 0.0, 1e-4, 0.0, 1e-4, 1e-4, 0.01]mol

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

# function update_vars(m, tstop)
    # m = @set m.records[1].vars = build_vars(Vars(), 1:ustrip(tstop))
    # m = @set m.records[2].vars = build_vars(Vars(), 1:ustrip(tstop))
    # m
# end

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
