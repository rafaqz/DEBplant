models[:maturity] = Plant(
    environment = tas,
    params = (
        Params(
            rate_formula = FZeroRate(),
            assimilation_pars = KooijmanWaterPotentialPhotosynthesis(
                potential_modifier = ZhouPotentialDependence(
                    s = 2.0u"MPa^-1",
                    ψ = -1.0u"MPa",
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
                M_Vref = 0.1u"mol",
                M_Vscaling = 1.0u"mol",
            ),
            allometry_pars = Allometry(
                β0 = 0.0024000000000000002u"g",
                β1 = 0.1u"m",
                α = 0.1,
            ),
            maturity_pars = Maturity(
                n_N_M = 0.1,
                j_E_mat_mai = 0.001u"d^-1",
                κmat = 0.05,
                threshold = 10.0u"mol",
            ),
            trans_pars = nothing,
            rejection_pars = LosslessRejection(),
            germination_pars = nothing,
            production_pars = nothing,
        ),
        Params(
            rate_formula = FZeroRate(),
            assimilation_pars = ConstantNAssim(
                uptake = 0.1u"μmol*mol^-1*s^-1",
            ),
            shape_pars = Plantmorph(
                M_Vref = 0.1u"mol",
                M_Vscaling = 1.0u"mol",
            ),
            allometry_pars = Allometry(
                β0 = 0.0024000000000000002u"g",
                β1 = 0.1u"m",
                α = 0.1,
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
            y_V_E = 0.7,
            y_E_EC = 0.7,
            y_E_EN = 30.0,
            n_N_V = 0.03,
            n_N_E = 0.025,
            w_V = 25.0u"g*mol^-1",
            w_N = 25.0u"g*mol^-1",
        ),
        feedback_pars = DissipativeAutophagy(
            r_EC_V = 0.0,
            r_EN_V = 0.5,
            K_autophagy = 1.0e-6u"mol",
        ),
        tempcorr_pars = TempCorrLowerUpper(
            reftemp = 300.0u"K",
            arrtemp = 2000.0u"K",
            tbelow = -30.0u"K",
            arrlower = 20000.0u"K",
            tabove = 5.0u"K",
            arrupper = 70000.0u"K",
        ),
        catabolism_pars = CatabolismCN(
            k_E = 0.2u"d^-1",
        ),
        maintenance_pars = Maintenance(
            j_E_mai = 0.01u"d^-1",
        ),
    ),
)