models[:bb] = Plant(
    environment = first(values(environments)),
    time = 0hr:1hr:8760hr*2,
    params = (
        Params(
            rate_formula = FZeroRate(),
            assimilation_pars = BBPotentialCAssim(
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
                    photosynthesis = FvCBPhotosynthesis(
                        flux = PotentialModifiedFlux(
                            flux = Flux(
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
                            potential_model = ZhouPotentialDependence(
                                s = 2.836u"MPa^-1",
                                ψ = -0.958u"MPa",
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
                            dayresp = 0.8,
                            rd0 = 0.01u"μmol*m^-2*s^-1",
                            tbelow = 173.15u"K",
                            tref = 298.15u"K",
                        ),
                        stomatal_conductance = BallBerryStomatalConductance(
                            g0 = 0.5u"μmol*m^-2*s^-1",
                            gs_submodel = BallBerryGSsubModel(
                                gamma = 0.0u"μmol*mol^-1",
                                g1 = 7.0,
                            ),
                            soilmethod = PotentialSoilMethod(
                                swpexp = 0.5u"kPa^-1",
                            ),
                        ),
                    ),
                ),
                SLA = 24.0u"m^2*kg^-1",
            ),
            shape_pars = Plantmorph(
                M_Vref = 0.001u"mol",
                M_Vscaling = 7.4u"mol",
            ),
            allometry_pars = Allometry(
                β1 = 0.09u"m",
                α = 0.2,
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
            shape_pars = Plantmorph(
                M_Vref = 0.0007u"mol",
                M_Vscaling = 1.0u"mol",
            ),
            allometry_pars = Allometry(
                β1 = 1.0u"m",
                α = 0.2,
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
            y_V_E = 1.0,
            y_E_EC = 0.9,
            y_E_EN = 30.0,
            n_N_V = 0.03,
            n_N_E = 0.025,
            w_V = 25.0u"g*mol^-1",
        ),
        resorption_pars = StructuralLossResorption(
            K_resorption = 1.0e-6,
        ),
        tempcorr_pars = ParentTardieu(
            ΔH_A = 63.5u"kJ*mol^-1",
            α = 3.5,
            t0 = 300.0u"K",
        ),
        catabolism_pars = CatabolismCNshared(
            k = 0.6u"d^-1",
        ),
        maintenance_pars = Maintenance(
            j_E_mai = 0.01u"d^-1",
        ),
    ),
)
