models[:potato] = Plant(
    environment = first(values(environments)),
    time = 0hr:1hr:8760hr*2,
    params = (
        Params(
            rate_formula = FZeroRate(),
            assimilation_pars = KooijmanWaterPotentialPhotosynthesis(
                potential_modifier = ZhouPotentialDependence(
                    s = 2.0u"MPa^-1",
                    ψ = -1.0u"MPa",
                ),
                vars = CarbonVars(
                ),
                k_C_binding = 10000.0u"μmol*mol^-1*s^-1",
                k_O_binding = 10000.0u"μmol*mol^-1*s^-1",
                K_C = 2.232142857142857e-6u"mol*L^-1",
                K_O = 9.375e-5u"mol*L^-1",
                J_L_K = 2000.0u"mol*m^-2*s^-1",
                j_L_Amax = 100.01u"μmol*m^-2*s^-1",
                j_C_Amax = 20.0u"μmol*m^-2*s^-1",
                j_O_Amax = 0.1u"μmol*m^-2*s^-1",
                SLA = 24.0u"m^2*kg^-1",
            ),
            shape_pars = Plantmorph(
                M_Vref = 0.2834948325853611u"mol",
                M_Vscaling = 33.36201074400118u"mol",
            ),
            allometry_pars = Allometry(
                β1 = 0.1u"m",
                α = 0.23101297000831605,
            ),
            maturity_pars = nothing,
            trans_pars = nothing,
            rejection_pars = LosslessRejection(),
            germination_pars = nothing,
            production_pars = nothing,
        ),
        Params(
            rate_formula = FZeroRate(),
            assimilation_pars = ConstantNAssim(
                uptake = 0.2u"μmol*mol^-1*s^-1",
            ),
            shape_pars = Plantmorph(
                M_Vref = 0.32595016692412887u"mol",
                M_Vscaling = 9.283177667225555u"mol",
            ),
            allometry_pars = Allometry(
                β1 = 1.321941148466029u"m",
                α = 0.13219411484660287,
            ),
            maturity_pars = nothing,
            trans_pars = nothing,
            rejection_pars = LosslessRejection(),
            germination_pars = nothing,
            production_pars = nothing,
        ),
    ),
    shared = SharedParams(
        su_pars = ParallelComplementarySU(),
        core_pars = DEBCore(
            y_V_E = 0.8,
            y_E_EC = 0.65000035,
            y_E_EN = 30.0,
            n_N_V = 0.03,
            n_N_E = 0.025,
            w_V = 25.0u"g*mol^-1",
        ),
        feedback_pars = StructuralLossAutophagy(
            K_autophagy = 8.697490026177835e-6u"mol",
        ),
        tempcorr_pars = ParentTardieu(
            ΔH_A = 63.5u"kJ*mol^-1",
            α = 3.5,
            t0 = 300.0u"K",
        ),
        catabolism_pars = CatabolismCNshared(
            k = 0.35u"d^-1",
        ),
        maintenance_pars = Maintenance(
            j_E_mai = 0.0049770235643321085u"d^-1",
        ),
    ),
)