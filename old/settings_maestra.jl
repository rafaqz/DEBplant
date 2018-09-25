using DynamicEnergyBudgetsBase
using DynamicEnergyBudgets
using DynamicEnergyBudgets: apply_maestra!
using PlantPhysiology
using NicheMap
using Unitful

macro inserttype(newtype, typenames...)
    fields = [:($(symbol("I_$i"))::T) for i=1:N]
    quote
        struct $(typename){T}
            $(fields...)
        end
    end
end

mutable struct MaestraPars{GS,CP,SD}
    k_E::typeof(0.2u"mol*mol^-1*d^-1")
    k_EC::typeof(0.2u"mol*mol^-1*d^-1")
    k_EN::typeof(0.2u"mol*mol^-1*d^-1")
    j_E_mai::typeof(0.001u"mol*mol^-1*d^-1")
    j_E_rep_mai::typeof(0.001u"mol*mol^-1*d^-1")
    j_P_mai::typeof(0.01u"mol*mol^-1*d^-1")
    SLA::typeof(9.10u"m^2*g^-1")
    K_autophagy::typeof(0.000001u"mol")
    K_C::typeof(40.0u"mol*L^-1")
    K_O::typeof(0.0021u"mol*L^-1")
    X_C::typeof(400.0u"mol*L^-1")
    X_O::typeof(0.21u"mol*L^-1")
    M_Vgerm::typeof(0.5u"mol")
    M_Vrep::typeof(10.0u"mol")
    M_Vref::typeof(4.0u"mol")
    M_Vscaling::typeof(400.0u"mol")
    κEC::typeof(0.2)
    κEN::typeof(0.5)
    κsoma::typeof(0.6)
    κrep::typeof(0.05)
    y_V_E::typeof(0.7u"mol*mol^-1")
    y_P_V::typeof(0.02u"mol*mol^-1")
    y_E_CH_NO::typeof(1.5u"mol*mol^-1")
    y_E_EN::typeof(0.5u"mol*mol^-1")
    y_E_ET::typeof(0.8u"mol*mol^-1")
    y_EN_ENT::typeof(1.0u"mol*mol^-1")
    n_N_P::typeof(0.0u"mol*mol^-1")
    n_N_V::typeof(0.15u"mol*mol^-1")
    n_N_EC::typeof(0.0u"mol*mol^-1")
    n_N_EN::typeof(10.0u"mol*mol^-1")
    n_N_E::typeof(0.2u"mol*mol^-1")
    w_P::typeof(25.0u"g*mol^-1")
    w_V::typeof(25.0u"g*mol^-1")
    w_EC::typeof(25.0u"g*mol^-1")
    w_EN::typeof(25.0u"g*mol^-1")
    w_E::typeof(25.0u"g*mol^-1")
    REFERENCE_TEMP::typeof(310.0u"K")
    ARRH_TEMP::typeof(2000.0u"K")
    LOWER_BOUNDARY::typeof(280.0u"K")
    ARRH_LOWER::typeof(20000.0u"K")
    UPPER_BOUNDARY::typeof(315.0u"K")
    ARRH_UPPER::typeof(70000.0u"K")
    j_NH_Amax::typeof(50.0u"μmol*s^-1")
    j_NO_Amax::typeof(50.0u"μmol*s^-1")
    K_NH::typeof(10.0u"mmol*L^-1")
    K_NO::typeof(10.0u"mmol*L^-1")
    K_H::typeof(10.0u"mol*L^-1")
    X_NH::typeof(5.0u"mmol*L^-1")
    X_NO::typeof(10.0u"mmol*L^-1")
    X_H::typeof(10.0u"mol*L^-1")
    ρNO::typeof(0.7)
    y_E_CH_NH::typeof(1.25u"mol*mol^-1")
    y_EC_ECT::typeof(1.0u"mol*mol^-1")
    photoparams::typeof(1.0u"mol*mol^-1")


    modelgs::GS
    compensation_point::CP
    soildata::SD
    itermax::typeof(100)
    rdfipt::typeof(1.0)
    tuipt::typeof(1.0)
    tdipt::typeof(1.0)
    rnet::typeof(0.0u"W*m^-2")
    windspeed::typeof(0.0u"m*s^-1")
    par::typeof(0.0u"mol*m^-2*s^-1")
    tair::typeof(0.0u"°C")
    pressure::typeof(101.0u"Pa")
    soilmoist::typeof(0.5)
    soilpot::typeof(0.5)
    soilpotshade::typeof(0.5)
    ca::typeof(400.0u"μmol*mol^-1")
    rh::typeof(1.0)
    vpd::typeof(1.5u"Pa")
    vpdmin::typeof(0.5u"Pa")
    jmax25::typeof(184.0u"μmol*m^-2*s^-1")
    vcmax25::typeof(110.0u"μmol*m^-2*s^-1")
    delsj::typeof(631.88u"J*mol^-1*K^-1")
    delsc::typeof(629.26u"J*mol^-1*K^-1")
    eavj::typeof(29680.0u"J*mol^-1")
    edvj::typeof(200000.0u"J*mol^-1")
    eavc::typeof(58550.0u"J*mol^-1")
    edvc::typeof(200000.0u"J*mol^-1")
    tvjup::typeof(0.0u"°C")
    tvjdn::typeof(0.0u"°C")
    theta::typeof(0.85)
    ajq::typeof(0.425)
    rd0::typeof(0.92u"μmol*m^-2*s^-1")
    q10f::typeof(1.92u"K^-1")
    rtemp::typeof(25.0u"°C")
    dayresp::typeof(1.0)
    tbelow::typeof(25.0u"°C")
    leafwidth::typeof(0.05u"m")
    gsc::typeof(1.0u"mol*m^-2*s^-1")
    smd1::typeof(1.0)
    smd2::typeof(1.0)
    swpexp::typeof(1.0)
    g0::typeof(0.0u"mol*m^-2*s^-1")
    gamma::typeof(0.0u"μmol*mol^-1")
    g1::typeof(4.0)
    # d0l::typeof(1.0u"u"Pa"") # Three-Par Ball-Berry only
    gk::typeof(0.5)
end


function build_settings(tspan; environment=[], use_environment=false, save_intermediate=false, timestep_days=1.0/24.0)
    state_type = StatePVMCNE
    if state_type == StatePVMCN
        u0 = [0.0, 1e-4, 0.0, 1e-4, 1e-4, 0.0, 1e-4, 0.0, 1e-4, 10.0]u"mol" # Initial value
        # u0 = [0.0, 1e-4, 1e-2, 1e-2, 0.0, 1e-4, 8.0, 2.0, 0.0, 1e-4, 8.0, 2.0] # Initial value
    elseif state_type == StatePVCNE
        u0 = [0.0, 1e-4, 1e-4, 1e-4, 1e-4, 0.0, 1e-4, 1e-4, 1e-4, 10.0]u"mol" # Initial value
        # u0 = [0.0, 1e-4, 1e-4, 1e-4, 1e-4, 0.0, 1e-4, 1e-4, 1e-4, 10.0, 0.0, 1e-4, 1e-4, 1e-4, 10.0] # Initial value
    elseif state_type == StatePVMCNE
        u0 = [0.0, 1e-4, 0.0, 1e-4, 1e-4, 1e-4, 0.0, 1e-4, 0.0, 1e-4, 1e-4, 10.0]u"mol" # Initial value
        # u0 = [0.0, 1e-4, 1e-4, 1e-4, 1e-4, 0.0, 1e-4, 1e-4, 1e-4, 10.0, 0.0, 1e-4, 1e-4, 1e-4, 10.0] # Initial value
    end

    load_settings(DEBSettings, MaestraPars, OrderedDict(
        :settings => OrderedDict(
            :u0 => u0,
            :tspan => tspan,
            :environment => environment,
            :use_environment => use_environment,
            :apply_environment! => apply_maestra!,
            :save_intermediate => save_intermediate,
            :timestep_days => timestep_days,
            :state_type => state_type
        ),
        :structures => OrderedDict(
            :Leaf => OrderedDict(
                :functions => OrderedDict(
                    :area => area_mass_kooijman,
                    :assim => shoot_assimilation!,
                    :assim_sub => maestra,
                    :rate => find_rate,
                ),
                :params => OrderedDict(
                    # J_L_F => (watts_to_light_mol(800.0) # mol/m²s, flux of useful photons
                    # TODO: should this acclimatise over time? this would require another state variable.

                    # J_L_F should feed back to affect plant state: photosynthesis as a sensory
                    # as well as energetic process [_@huner1998energy
                    # J_L_K => (watts_to_light_mol(300.0), "mol/m²s, half-saturation flux of useful photons"
                    # Max specific uptake parameters that relate uptake to active surface area
                    # These are usually measured in μmol m⁻² s⁻¹ (wikipedia:Photosynthetic capacity).
                    # From [_@walker2014relationship :...full range of photosynthetically active radiation
                    # (PAR 0–1500 μmol·m−2·s−1) and three levels of Vcmax (25, 50 & 90 μmol·m−2·s−1)
                    # j_L_Amax => (μmol_to_mol(20.0) #[@2006seasonality] umol.m⁻².s⁻¹, max spec uptake of useful photons
                    # j_C_Amax => (μmol_to_mol(90.0), "mol.m⁻².s⁻¹, max spec uptake of carbon dioxide"
                    # j_O_Amax => (μmol_to_mol(0.001), "mol.m⁻².s⁻¹, max spec uptake of oxygen"
                    # Binding rates of gases to quantify photo-respiration
                    # Is this braodly similar accross plants?
                    # k_C_binding => (1.0, :temp, "mols.s⁻¹, scaling rate for carbon dioxide "
                    # k_O_binding => (1.0, :temp, "mols.s⁻¹, scaling rate for oxygen"
                    # Turnover. This controls metabolic rate along with area/mass relation and
                    # current reserves.
                    :k_E => (0.2u"mol*mol^-1*d^-1", (0.0u"mol*mol^-1*d^-1", 1.0u"mol*mol^-1*d^-1"), (:exposed, :time, :temp), "shoots reserve turnover rate"),
                    :k_EC => (0.2u"mol*mol^-1*d^-1", (0.0u"mol*mol^-1*d^-1", 1.0u"mol*mol^-1*d^-1"), (:exposed, :time, :temp), "shoots C-reserve turnover rate"),
                    :k_EN => (0.2u"mol*mol^-1*d^-1", (0.0u"mol*mol^-1*d^-1", 1.0u"mol*mol^-1*d^-1"), (:exposed, :time, :temp), "shoots N-reserve turnover rate"),
                    # Specific maintenancs costs in terms of reserve fluxes.
                    # Must be per day not s, or these numbers are crazy.
                    # - these costs are paid to maintain structural mass and reproductive maturity
                    :j_E_mai => (0.001u"mol*mol^-1*d^-1", (0.0u"mol*mol^-1*d^-1", 0.01u"mol*mol^-1*d^-1"), (:exposed, :time, :temp), "shoots spec somatic maint costs."),
                    :j_E_rep_mai => (0.001u"mol*mol^-1*d^-1", (0.0u"mol*mol^-1*d^-1", 0.01u"mol*mol^-1*d^-1"), (:exposed, :time, :temp), "shoots spec maturity maint costs "),
                    # Production parameters
                    # - these parameters play no dynamic role but can dominate weights
                    :j_P_mai => (0.01u"mol*mol^-1*d^-1", (0.0u"mol*mol^-1*d^-1", 0.1u"mol*mol^-1*d^-1"), (:exposed, :time, :temp), "shoot product formation linked to maintenance"),

                    :SLA => (9.10u"m^2*g^-1", (5.0u"m^2*g^-1", 30.0u"m^2*g^-1"), "Ferns 17.4 Forbs 26.2 Graminoids 24.0 Shrubs 9.10 Trees 8.30"),

                    :K_autophagy => (0.000001u"mol", (0.0000001u"mol", 0.00001u"mol")),
                    :K_C => (fraction_per_litre_gas_to_mols(40.0/1e6)u"mol*L^-1", "half-saturation concentration of carbon dioxide"),
                    :K_O => (fraction_per_litre_gas_to_mols(0.0021)u"mol*L^-1", "half-saturation concentration of oxygen"),
                    :X_C => (fraction_per_litre_gas_to_mols(400.0/1e6)u"mol*L^-1", "carbon dioxide @ 400ppm"),
                    :X_O => (fraction_per_litre_gas_to_mols(0.21)u"mol*L^-1", "oxygen (21% volume in air) "),
                    # Life stage parameters
                    :M_Vgerm => (0.5u"mol", "shoots structural mass at germination"),
                    :M_Vrep => (10.0u"mol", "shoots structural mass at start reproduction"), # TODO: isn't this variable/seasonally triggered?

                    # Parameters that link active surface area to structural mass
                    # - they describe the development through V1- iso- and V0-morphs
                    # TODO: estimate these from functional traints SLA/LMA and/or
                    # stem/branch wood specific gravity and final height?
                    # The curve function probably also needs replacing.
                    :M_Vref => (4.0u"mol", (0.4u"mol", 20.0u"mol"), (:exposed,)),
                    :M_Vscaling => (400.0u"mol", (40.0u"mol", 2000.0u"mol"), (:exposed,), "shoots scaling mass"),

                    # Partitioning parameters: dimensionless fractions
                    :κEC => (0.2, (0.0, 1.0), (:exposed,), "shoots non-processed C-reserve returned to C-reserve,"),
                    #    the remaining fraction is translocated to the root
                    :κEN => (0.5, (0.0, 1.0), (:exposed), "shoots non-processed N-reserve returned to N-reserve"),
                    :κsoma => (0.6, (0.0, 1.0), (:exposed,), "shoots reserve flux allocated to soma"),
                    :κrep => (0.05, (0.0, 1.0), (:exposed,), "shoots reserve flux allocated to development/reprod."),
                    :y_V_E => (0.7u"mol*mol^-1", (:exposed), "from shoots reserve to structure: 0.3 lost as CO2?"),
                    :y_P_V => (0.02u"mol*mol^-1", (:exposed,), "shoot product formation linked to growth: where does this carbon come from?"),
                    :y_E_CH_NO => (1.5u"mol*mol^-1", (:exposed,), "from shoots C-reserve to reserve, using nitrate: 0.75 EC per E."),
                    :y_E_EN => (0.5u"mol*mol^-1", (:exposed,), "from shoots N-reserve to reserve: 2 EN per E"),
                    :y_E_ET => (0.8u"mol*mol^-1", (:exposed,), "from shoots reserve to roots reserve:"),
                    :y_EN_ENT => (1.0u"mol*mol^-1", (:exposed,), "from roots N-reserve to shoots N-reserve"),

                    # Chemical indices: elemental_nitrogen/elemental_carbon
                    # - only n_N_EN and n_N_E play a dynamic role
                    # - the others are only used to evaluate the nitrogen balance
                    :n_N_P => (0.0u"mol*mol^-1", "N/C in product (wood)"),
                    :n_N_V => (0.15u"mol*mol^-1", "N/C in structure. Shouldnt this be identical to the reserve?"),
                    :n_N_M => (10.0u"mol*mol^-1", "N/C in M-reserve"),
                    :n_N_EC => (0.0u"mol*mol^-1", "N/C in C-reserve"),
                    :n_N_EN => (10.0u"mol*mol^-1", "N/C in N-reserve"),
                    :n_N_E => (0.2u"mol*mol^-1", "N/C in reserve. This should be calculated, not constant. 1.8181??? (10/11 * 0.3)/1.5"),

                    # Parameters that link moles to grams (wet weight)
                    # These should be calculated from first principles, not preset.
                    :w_P => (25.0u"g*mol^-1", "mol-weight of shoot product (wood)"),
                    :w_V => (25.0u"g*mol^-1", "mol-weight of shoot structure"),
                    :w_M => (25.0u"g*mol^-1", "mol-weight of shoot matureity reserve:"),
                    :w_EC => (25.0u"g*mol^-1", "mol-weight of shoot C-reserve"),
                    :w_EN => (25.0u"g*mol^-1", "mol-weight of shoot N-reserve"),
                    :w_E => (25.0u"g*mol^-1", "mol-weight of shoot reserve"),

                    :REFERENCE_TEMP => (310.0u"K", (273.0u"K", 325.0u"K"), "temp for which rate pars are given"),
                    :ARRH_TEMP => (2000.0u"K", (200.0u"K", 4000.0u"K"), "Arrhenius temp"),
                    :LOWER_BOUNDARY => (280.0u"K", (273.0u"K", 325.0u"K"), "lower boundary tolerance range"),
                    :ARRH_LOWER => (20000.0u"K", (2000.0u"K", 40000.0u"K"), "Arrhenius temp for lower boundary"),
                    :UPPER_BOUNDARY => (315.0u"K", (273.0u"K", 325.0u"K"), "upper boundary tolerance range"),
                    :ARRH_UPPER => (70000.0u"K", (7000.0u"K", 140000.0u"K"), "Arrhenius temp for upper boundary"),

                    # :LEAF_D => (0.01, (0.0001, 0.1)),
                    # :LEAF_W => (0.01, (0.0001, 0.1)),
                    # :LEAF_EMISSIVITY => (0.5, (0.0, 1.0)),

                    # Vpm25 => (-1.0
                    # TPU25 => (15.0
                    # Rd25 => (1.0
                    # θ => (0.7
                    # EaVc => (64800.0
                    # Eaj => (37000.0
                    # Hj => (220000.0
                    # Sj => (710.0
                    # Hv => (219400.0
                    # EaVp => (-1.0
                    # Sv => (-1.0
                    # Eap => (47100.0
                    # Ear => (66400.0
                    # g0 => (0.02
                    # g1 => (10.0
                    # stoma_ratio => (0.5
                    # leaf_width => (0.05
                    # leaf_angfact => (1.0
                    # photo_flux_density => (1700.0
                    # Tair => (25.0
                    # CO2 => (370.0
                    # rel_humidity => (65.0
                    # wind => (0.8
                    # pressure => (1.0

                    :modelgs => (Medlyn()),
                    :compensation_point => (Bernacchi()),
                    :soildata => (Potential()),
                    :itermax => (100),
                    :rdfipt => (1.0),
                    :tuipt => (1.0),
                    :tdipt => (1.0),
                    :rnet => (0.0u"W*m^-2"),
                    :windspeed => (0.0u"m*s^-1"),
                    :par => (0.0u"mol*m^-2*s^-1"),
                    :tair => (0.0u"°C"),
                    :pressure => (101.0u"Pa"),
                    :soilmoist => (0.5),
                    :ca => (400.0u"μmol*mol^-1"),
                    :rh => (1.0), # Fraction. Only for Ball-Berry. This value is made up.
                    :vpd => (1.5u"Pa"),
                    :vpdmin => (0.5u"Pa"),
                    :jmax25 => (184.0u"μmol*m^-2*s^-1"),
                    :vcmax25 => (110.0u"μmol*m^-2*s^-1"),
                    :delsj => (631.88u"J*mol^-1*K^-1"),
                    :delsc => (629.26u"J*mol^-1*K^-1"), # different to [@black2009Carbon] - needs to be per K!
                    :eavj => (29680.0u"J*mol^-1"),
                    :edvj => (200000.0u"J*mol^-1"),
                    :eavc => (58550.0u"J*mol^-1"),
                    :edvc => (200000.0u"J*mol^-1"),
                    :tvjup => (1.0u"°C"),
                    :tvjdn => (1.0u"°C"),
                    :theta => (0.85),
                    :ajq => (0.4), # quantum yield of electron transport
                    :rd0 => (0.92u"μmol*m^-2*s^-1"),
                    :q10f => (1.92u"K^-1"),
                    :rtemp => (25.0u"°C"),
                    :dayresp => (1.0), # fraction made up
                    :tbelow => (25.0u"°C"), # made up
                    :leafwidth => (0.05u"m"),
                    :gsc => (1.0u"mol*m^-2*s^-1"),
                    :smd1 => (1.0), # for Soildata =>=> Deficit
                    :smd2 => (1.0),
                    :swpexp => (1.0), # for Soildata =>=> Potential
                    :g0 => (0.1u"mol*m^-2*s^-1"), # All Ball-Berry
                    :gamma => (0.1u"μmol*mol^-1"),
                    :g1 => (4.0), # dimensionless # Ball-Berry-Leuning only
                    # d0l::typeof(1.0u"Pa") # Three-Par Ball-Berry only
                    :gk => (0.5), # exponent
                )
            ),
            :Root => OrderedDict(
                :functions => OrderedDict(
                    :area => area_mass_kooijman,
                    :assim => root_assimilation!,
                    :assim_sub => nitrogen_uptake,
                    :rate => find_rate,
                ),
                :params => OrderedDict(
                    :j_NH_Amax => (50.0u"μmol*s^-1", (0.1u"μmol*s^-1", 1000.0u"μmol*s^-1"), "max spec uptake of ammonia"), 
                    :j_NO_Amax => (50.0u"μmol*s^-1", (0.1u"μmol*s^-1", 1000.0u"μmol*s^-1"), "max spec uptake of nitrate"), 
                    :j_E_rep_mai => (0.0u"mol*mol^-1*d^-1", "roots spec maturity maint costs"),

                    # Water affects root nutrient uptake via the saturation parameters
                    :K_NH => (10.0u"mmol*L^-1", "half-saturation concentration of ammonia"),
                    :K_NO => (10.0u"mmol*L^-1", "half-saturation concentration of nitrate"),
                    :K_H => (10.0u"mol*L^-1", "half-saturation concentration of water"),
                    :X_NH => (5.0u"mmol*L^-1", "ammonia"),
                    :X_NO => (10.0u"mmol*L^-1", "concentration of nitrate see e.g. [_@crawford1998molecular]"),
                    :X_H => (10.0u"mol*L^-1"),

                    # Life stage parameters
                    :M_Vgerm => (0.3u"mol", "roots structural mass at germination"),
                    :M_Vrep => (10.0u"mol", "roots structural mass at start reproduction"),

                    :κEC => (0.5, (0.0, 1.0), (:exposed,), "roots  non-processed C-reserve returned to C-reserve"),
                    #    the remaining fraction is translocated to the shoot
                    :κEN => (0.2, (0.0, 1.0), (:exposed,), "roots  non-processed N-reserve returned to N-reserve"),
                    :κsoma => (0.5, (0.0, 1.0), (:exposed,), "roots  reserve flux allocated to soma"),
                    :κrep => (0.0, (0.0, 1.0), "shoots reserve flux allocated to development/reprod."),

                    :ρNO => (0.7, (0.0, 1.0), "weights preference for nitrate relative to ammonia."), # 1 or less but why?

                    # Yield coefficients (see also production parameters)
                    :y_E_CH_NO => (1.5u"mol*mol^-1", "from roots C-reserve to reserve, using nitrate"),
                    :y_E_CH_NH => (1.25u"mol*mol^-1", "from roots C-reserve to reserve, using ammonia"),
                    :y_EC_ECT => (1.0u"mol*mol^-1", "from shoots C-reserve to roots C-reserve"),
                    :y_E_EN => (0.3u"mol*mol^-1", "from roots N-reserve to reserve: why is this different to shoots?"),
                )
            )
        )
   ))
end
