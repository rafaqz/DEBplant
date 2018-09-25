function catabolism!(o, u::AbstractStateE, t::Number)
    p = o.params; v = o.vars; J1 = o.J1
    scaledturnover = (p.k_E,) .* v.scale
    m = u.E / u.V
    v.rate = find_rate(v, (m, scaledturnover, p.j_E_mai, p.y_V_E, p.κsoma))
    J1[:E,:cat] = catabolic_fluxes(ureserve, scaledturnover, v.rate)
    return nothing
end
function catabolism!(o, u::AbstractStateCN, t::Number)
    p = o.params; v = o.vars; J1 = o.J1
    scaledturnover = (p.k_EC, p.k_EN) .* v.scale
    ureserve = (u.C, u.N)
    m = ureserve ./ u.V
    v.rate = find_rate(v, (m, scaledturnover, p.j_E_mai, p.y_E_CH_NO, p.y_E_EN, p.y_V_E, p.κsoma))
    (J1[:C,:cat], J1[:N,:cat]) = catabolic_fluxes(ureserve, scaledturnover, v.rate)
    (J1[:C,:rej], J1[:N,:rej], J1[:E,:cat]) =
        synthesizing_unit(J1[:C,:cat], J1[:N,:cat], p.y_E_CH_NO, p.y_E_EN)
    return nothing
end

function translocate!(o, on, u::AbstractStateE)
    κtra = calc_κtra(o.params)
    κtraT = calc_κtra(on.params)
    # Translocation of C, N and general reserves between organs
    Cgain = κtraT * on.state.J1[:E,:cat]
    drain = -κtra * o.J1[:E,:cat]
    o.J[:E,:tra] += Cdrain + Cgain
    o.J1[:E,:los] -= gain
    return nothing
end
function translocate!(o, on, u::AbstractStateCN)
    J = o.J; J1 = o.J1;
    Jn = on.state.J; J1n = on.state.J1;
    κtra = calc_κtra(o.params)
    κtraT = calc_κtra(on.params)

    # Translocation of C, N and general reserves between organs
    tra_C = -κtra * J1[:C,:cat]
    tra_N = -κtra * J1[:N,:cat]

    J[:C,:tra] = tra_C * 1/o.paramo.y_E_CH_NO
    Cdrain = -κtra * J1[:C,:cat]
    Cgain = o.params.y_E_ET * κtraT * J1n[:C,:cat]
    J[:C,:tra] = Cdrain + Cgain

    J[:N,:tra] = tra_N * 1/o.params.y_E_EN
    Ndrain = -κtra * J1[:N,:cat]
    Ngain = o.params.y_E_ET * κtraT * J1n[:N,:cat]
    J[:N,:tra] = Ndrain + Ngain
    return nothing
end

function assimilation!(o, on, f::CarbonAssimilation, u::AbstractStateCN)::Nothing
    p = o.params; v = o.vars
    if !germinated(u.V, p.M_Vgerm) 
        on.J[:N,:rej] = o.J[:C,:ass] = o.J[:N,:ass] = zero(o.J[1,1])
        return nothing
    end
    J1_EC_ass = photosynthesize(f, o) * u.V * v.scale

    # Rejected N-reserve from root
    J1_EN_ass = -p.y_EN_ENT * on.J1[:N,:rej] 

    (o.J[:C,:ass], o.J[:N,:ass]) = (J1_EC_ass, J1_EN_ass)
    return nothing
end


# function assimilation!(o, on, f::NH4_NO3_Assimilation, u::AbstractStateCN)::Nothing
#     p = o.params; v = o.vars
#     if !germinated(u.V, p.M_Vgerm) 
#         on.J[:C,:rej] = o.J[:C,:ass] = o.J[:N,:ass] = zero(o.J[1,1]) 
#         return nothing
#     end
#     # Arriving nitrogen
#     (J_N_ass, J_NO_ass, J_NH_ass) = uptake_nitrogen(f.formulation, o, on)

#     # Rejected C-reserve from shoot
#     J1_EC_ass = -p.y_EC_ECT * on.J1[:C,:rej]

#     θNH = J_NH_ass/J_N_ass                           # Fraction of ammonia in arriving N-flux
#     θNO = 1 - θNH                                    # Fraction of nitrate in arriving N-flux
#     y_E_CH = θNH * f.y_E_CH_NH + θNO * p.y_E_CH_NO   # Yield coefficient from C-reserve to reserve

#     (o.J[:C,:ass], o.J[:N,:ass]) = (J1_EC_ass, J1_EN_ass)

#     # Unused NH₄ remainder is lost so we recalculate N assimilation for NO₃ only
#     # How to deal with this without general reserve?
#     # o.J[:N,:ass] = (J_NO_ass - θNO * p.n_N_E * o.J[:E,:ass]) * 1/p.n_N_EN
#     return nothing
# end

# function uptake_nitrogen(f::KooijmanSLA_NH4_NO3_Assimilation, o, on)
#     p = o.params; v = o.vars
#     SAX_H = f.X_H * on.vars.scale # mol.L⁻¹ Why A_S? Evapotranspiration based on shoot area, not root area. See p. 199 of DEB book.
#     RAK_H = f.K_H * v.scale # mol.L⁻¹ TODO: only partial root area is not involved in N uptake [@robinson1991what].
#     K1_NH = half_saturation(RAK_H, SAX_H, f.K_NH) # mol.L⁻¹ Ammonia saturation
#     K1_NO = half_saturation(RAK_H, SAX_H, f.K_NO) # mol.L⁻¹ Nitrate saturation
#     j1_NH = half_saturation(f.X_NH, K1_NH, f.j_NH_Amax) # Arriving ammonia mols.m⁻².s⁻¹ 
#     j1_NO = half_saturation(f.X_NO, K1_NO, f.j_NO_Amax) # Arriving nitrate mols.m⁻².s⁻¹
#     area = o.state.V * v.scale * p.w_V * f.SLA
#     J1_NH_ass = j1_NH / area
#     J_NO_ass = j1_NO / area
#     J_N_ass = J1_NH_ass + f.ρNO * J_NO_ass # Total arriving N flux
#     return (J_N_ass, J_NO_ass, J1_NH_ass)
# end

feedback!(o, f::Autophagy, u::AbstractStateE) = begin
    hs = half_saturation(oneunit(f.K_autophagy), f.K_autophagy, o.vars.rate)
    autophagy = u.V * (oneunit(hs) - hs)
    o.J[:E,:gro] += autophagy
    o.J[:V,:gro] -= autophagy
    nothing
end

function reserve_loss!(J1, J, u::AbstractStateE, col, θloss)
    J1[:E,:los] -= J[:E,col] * θloss
    return nothing
end
function reserve_loss!(J1, J, u::AbstractStateCN, col, θloss)
    J1[:C,:los] -= J[:C,col] * θloss
    J1[:N,:los] -= J[:N,col] * θloss
    return nothing
end

function reserve_drain!(J, u::AbstractState, col, drain, θE, params)
    return nothing
end
function reserve_drain!(J, u::AbstractStateE, col, drain, θE, params)
    J[:E,col] = drain * θE
    return nothing
end
function reserve_drain!(J, u::AbstractStateCN, col, drain, θE, params)
    J[:C,col] = drain/params.y_E_CH_NO
    J[:N,col] = drain/params.y_E_EN
    return nothing
end

"""
Get bracket of possible rates with CN or CNE rserves
Calculated for y_E_CH_NO and y_E_EN -> ∞
"""
function rate_bracket(ureserve::NTuple{2}, A_turnover::NTuple{2}, args::Vararg)::NTuple{2}
    rate_bracket((ureserve[1], ureserve[2], 0.0), (A_turnover[1], A_turnover[2], 0.0), args...)
end
function rate_bracket(ureserve::NTuple{3}, A_turnover::NTuple{3}, 
                     j_E_mai, y_E_CH_NO, y_E_EN, y_V_E, κsoma)::NTuple{2}
    # Calculate the limits of the rate_formula function when C or N → ∞
    (uEC, uEN, uE) = ureserve
    (AEC, AEN, AE) = A_turnover
    lim(y) = (uE * AE + uEC * AEC - j_E_mai/(κsoma * y))/(1/(y_V_E * κsoma * y) + uE + uEC * y)
    (lim(y_E_EN), lim(y_E_CH_NO))
end
