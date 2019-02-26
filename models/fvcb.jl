models[:fvcb] = Plant(
    environment = first(values(environments)),
    params = (
        Params(
            rate_formula = FZeroRate(),
            assimilation_pars = FvCBPhotosynthesis(
                photoparams = FvCBEnergyBalance(
                    radiation_conductance = YingPingRadiationConductance(
                        rdfipt = 1.0,
                        tuipt = 1.0,
                        tdipt = 1.0,
                    ),
                    boundary_conductance = BoundaryConductance(
                        leafwidth = 0.05u"m",
                        gsc = 1.0u"mol*m^-2*s^-1",
                    ),
                    decoupling = McNaughtonJarvisDecoupling(),
                    photo = BallBerryModel(
                        gsmodel = BallBerryStomatalConductance(
                            gamma = 0.0u"μmol*mol^-1",
                            g1 = 7.0,
                        ),
                        soilmethod = ConstantSoilMethod(
                            soildata = NoSoilData(),
                        ),
                        g0 = 0.03u"mol*m^-2*s^-1",
                        vcjmax = VcJmax(
                            jmaxformulation = Jmax(
                                jmax25 = 184.0u"μmol*m^-2*s^-1",
                                delsj = 640.02u"J*K^-1*mol^-1",
                                eavj = 37259.0u"J*mol^-1",
                                edvj = 200000.0u"J*mol^-1",
                            ),
                            vcmaxformulation = NoOptimumVcmax(
                                vcmax25 = 110.0u"μmol*m^-2*s^-1",
                                eavc = 47590.0u"J*mol^-1",
                            ),
                        ),
                        compensation = BernacchiCompensation(
                            Kc25 = 404.9u"μmol*mol^-1",
                            Ko25 = 278400.0u"μmol*mol^-1",
                            Γ☆25 = 42.75u"μmol*mol^-1",
                            ΔHa_Kc = 79.43u"kJ*mol^-1",
                            ΔHa_Ko = 36.38u"kJ*mol^-1",
                            ΔHa_Γ☆ = 37.83u"kJ*mol^-1",
                            tref = 298.15u"K",
                        ),
                        rubisco_regen = RubiscoRegen(
                            theta = 0.4,
                            ajq = 0.324,
                        ),
                        respiration = Respiration(
                            q10f = 0.67u"K^-1",
                            dayresp = 1.0,
                            rd0 = 0.9u"μmol*m^-2*s^-1",
                            tbelow = 173.15u"K",
                            tref = 298.15u"K",
                        ),
                    ),
                ),
                SLA = 24.0u"m^2*kg^-1",
            ),
            shape_pars = Plantmorph(
                M_Vref = 0.1u"mol",
                M_Vscaling = 32.59501669241289u"mol",
            ),
            allometry_pars = Allometry(
                β0 = 0.0024000000000000002u"g",
                β1 = 0.1u"m",
                α = 0.1,
            ),
            maturity_pars = Maturity(
                j_E_mat_mai = 0.013u"d^-1",
                κmat = 0.36,
                threshold = 0.3926081300080543u"mol",
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
                M_Vscaling = 37.47634845720768u"mol",
            ),
            allometry_pars = Allometry(
                β0 = 0.0024000000000000002u"g",
                β1 = 1.5199110829529336u"m",
                α = 0.21049041445120198,
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
            reftemp = 303.68u"K",
            arrtemp = 2000.0u"K",
            tbelow = -30.0u"K",
            arrlower = 20000.0u"K",
            tabove = 5.0u"K",
            arrupper = 70000.0u"K",
        ),
        catabolism_pars = CatabolismCN(
            k_E = 0.29u"d^-1",
        ),
        maintenance_pars = Maintenance(
            j_E_mai = 0.001788649529057435u"d^-1",
        ),
    ),
)