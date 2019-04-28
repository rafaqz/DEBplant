models[:emax] = Plant(
    environment = first(values(environments)),
    time = 0hr:1hr:8760hr*2,
    params = (
        Params(
            rate_formula = FZeroRate(),
            assimilation_pars = EmaxCAssim(
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
                    evapotranspiration = PenmanMonteithEvapotranspiration(),
                    photosynthesis = EmaxPhotosynthesis(
                        plantk = 3.0u"mmol*m^-2*MPa^-1*s^-1",
                        totsoilres = 0.5u"m^2*MPa*s*mmol^-1",
                        gsshape = HardMinimumGS(),
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
                        gsmodel = BallBerryStomatalConductance(
                            gamma = 0.0u"μmol*mol^-1",
                            g1 = 7.0,
                        ),
                        soilmethod = EmaxSoilMethod(
                            soilmethod = ConstantSoilMethod(
                                soildata = NoSoilData(),
                            ),
                            non_stomatal = ZhouPotentialDependence(
                                s = 2.0u"MPa^-1",
                                ψ = -1.0u"MPa",
                            ),
                        ),
                    ),
                ),
                SLA = 24.0u"m^2*kg^-1",
            ),
            shape_pars = Plantmorph(
                M_Vref = 0.00026438822969320576u"mol",
                M_Vscaling = 190.90969133236686u"mol",
            ),
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
                n_uptake = 0.2u"μmol*mol^-1*s^-1",
            ),
            shape_pars = Plantmorph(
                M_Vref = 0.001410960462143729u"mol",
                M_Vscaling = 53.121755658933736u"mol",
            ),
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
            y_V_E = 0.7,
            y_E_EC = 0.6900003100000001,
            y_E_EN = 30.0,
            n_N_V = 0.03,
            n_N_E = 0.025,
            w_V = 25.0u"g*mol^-1",
        ),
        feedback_pars = StructuralLossAutophagy(
            K_autophagy = 2.364489412645407e-6,
        ),
        tempcorr_pars = ParentTardieu(
            ΔH_A = 63.1u"kJ*mol^-1",
            α = 3.07,
            t0 = 297.96u"K",
        ),
        catabolism_pars = CatabolismCNshared(
            k = 0.7u"d^-1",
        ),
        maintenance_pars = Maintenance(
            j_E_mai = 0.010476157527896646u"d^-1",
        ),
    ),
)