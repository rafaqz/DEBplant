using Setfield
using Unitful, Plots, UnitfulRecipes, Microclimate, DataStructures, Flatten, FieldMetadata, 
      OrdinaryDiffEq, Photosynthesis, DynamicEnergyBudgets, Dates,
      DataStructures, ColorSchemes, JLD2, DimensionalData, FileIO

using DynamicEnergyBudgets: scaling_correction, define_organs,
      photosynthesis, split_state, HasCN, HasCNE, has_reserves, tempcorr_pars,
      assimilation_pars, parconv, w_V, build_vars, allometry_pars, tempcorr

using Unitful: °C, K, Pa, kPa, MPa, J, kJ, W, L, g, kg, g, mg, cm, m, s, hr, d, mol, mmol, μmol, σ

import Plots:px, pct, GridLayout

const STATELABELS = ["VS", "CS", "NS", "VR", "CR", "NR"]

function transect_from_saved(projectdir)
    locationspath = joinpath(projectdir, "data/locations.jld2")
    @load locationspath t1 t2 t3
    environments = OrderedDict{Symbol,Any}(:t1 => t1, :t2 => t2, :t3 => t3)
    env = t1
    tspan = (0:1:length(radiation(env)) - 1) * u"hr"
    environments, tspan
end

function transect_from_netcdf(microclimdir::String, years, shade)
    envgrid = load_grid(microclimdir, years, shade);

    t1 = MicroclimPoint(envgrid, CartesianIndex(65, 35))
    t2 = MicroclimPoint(envgrid, CartesianIndex(60, 35))
    r3 = MicroclimPoint(envgrid, CartesianIndex(55, 35))

    save(joinpath(dir, "data/locations.jld2"), Dict("t1" => t1, "t2" => t2, "t3" => t3))

    tspan = (0:1:length(radiation(env)) - 1) * u"hr"
    environments = OrderedDict{Symbol,Any}(:t1 => t1, :t2 => t2, :t3 => t3)
    environments, tspan
end

init_state(model::AbstractOrganism) = init_state(has_reserves.(define_organs(model, 1hr)), model)
init_state(::NTuple{2,HasCN}, model) = begin
    xdim = dims(first(model.records).J, X)
    A = zeros(length(xdim) * length(model.records)) * mol
    newxval = Val((:VS, :CS, :NS, :VR, :CR, :NR))
    newxdim = X(newxval, Categorical(), nothing)
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

