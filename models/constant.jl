models[:constant] = Plant(
    environment = environment,
    vars = vars,
    time = 0hr:1hr:8760hr*2,
    params = (
        Params(
            assimilation_pars = ConstantCAssim(
                c_uptake = 1.2u"μmol*mol^-1*s^-1",
            ),
            scaling_pars = Plantmorph(
                M_Vref = 0.02u"mol",
                M_Vscaling = 0.3u"mol",
            ),
            allometry_pars = Allometry(
                β1 = 0.2u"m",
                α = 0.2,
            ),
            maturity_pars = nothing,
            activetrans_pars = nothing,
            passivetrans_pars = LosslessPassiveTranslocation(),
            germination_pars = nothing,
            production_pars = nothing,
        ),
        Params(
            assimilation_pars = ConstantNAssim(
                n_uptake = 0.10u"μmol*mol^-1*s^-1",
            ),
            scaling_pars = Plantmorph(
                M_Vref = 0.02u"mol",
                M_Vscaling = 0.2u"mol",
            ),
            allometry_pars = Allometry(
                β1 = 1.0u"m",
                α = 0.2,
            ),
            maturity_pars = nothing,
            activetrans_pars = nothing,
            passivetrans_pars = LosslessPassiveTranslocation(),
            germination_pars = nothing,
            production_pars = nothing,
        ),
    ),
    shared = SharedParams(
        su_pars = ParallelComplementarySU(),
        core_pars = DEBCore(
            j_E_mai = 0.01u"d^-1",
            y_V_E = 1.0,
            y_E_EC = 0.9,
            y_E_EN = 30.0,
            n_N_V = 0.03,
            n_N_E = 0.025,
            w_V = 25.0u"g*mol^-1",
        ),
        resorption_pars = LosslessResorption(
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
    ),
)
