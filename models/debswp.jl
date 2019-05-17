models[:debswp] = Plant(
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
            shape_pars = Isomorph(),
            allometry_pars = Allometry(
                β1 = 0.093260334688322u"m",
                α = 0.19179102616724886,
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
                n_uptake = 0.15u"μmol*mol^-1*s^-1",
            ),
            shape_pars = Isomorph(),
            allometry_pars = Allometry(
                β1 = 1.0u"m",
                α = 0.19179102616724886,
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
            y_V_E = 0.5,
            y_E_EC = 0.51000049,
            y_E_EN = 30.0,
            n_N_V = 0.03,
            n_N_E = 0.025,
            w_V = 25.0u"g*mol^-1",
        ),
        resorption_pars = StructuralLossResorption(
            K_resorption = 1.047615752789665e-6,
        ),
        tempcorr_pars = ParentTardieu(
            ΔH_A = 63.5u"kJ*mol^-1",
            α = 3.5,
            t0 = 300.0u"K",
        ),
        catabolism_pars = CatabolismCNshared(
            k = 0.75u"d^-1",
        ),
        maintenance_pars = Maintenance(
            j_E_mai = 0.01519911082952934u"d^-1",
        ),
    ),
)
