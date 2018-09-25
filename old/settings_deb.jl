using DynamicEnergyBudgetsBase
using DynamicEnergyBudgets
using PlantPhysiology
using NicheMap
using NamedTuples
using Unitful
using UnitfulMoles
using Unitful: K, °C, W, J, m, s, Pa, mol, μmol
@mol Quanta 0.0

mutable struct DEBPars{GS,CP,SD}
    J_L_F::typeof(800.0u"mol*m^2*s^-1")            
    J_L_K::typeof(300.0u"mol*m^2*s^-1")
    j_L_Amax::typeof(20.0u"μmol*m^2*s^-1")
    j_C_Amax::typeof(90.0u"μmol*m^2*s^-1")
    j_O_Amax::typeof(0.001u"μmol*m^2*s^-1")
    k_C_binding::typeof(1.0u"mol*s^-1")
    k_O_binding::typeof(1.0u"mol*s^-1")
    k_E::typeof(0.2u"mol*d^-1")
    k_EC::typeof(0.2u"mol*d^-1")
    k_EN::typeof(0.2u"mol*d^-1")
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
    w_P::typeof(25.0u"mol*g^-1")
    w_V::typeof(25.0u"mol*g^-1")
    w_EC::typeof(25.0u"mol*g^-1")
    w_EN::typeof(25.0u"mol*g^-1")
    w_E::typeof(25.0u"mol*g^-1")
    REFERENCE_TEMP::typeof(310.0K)
    ARRH_TEMP::typeof(2000.0K)
    LOWER_BOUNDARY::typeof(280.0K)
    ARRH_LOWER::typeof(20000.0K)
    UPPER_BOUNDARY::typeof(315.0K)
    ARRH_UPPER::typeof(70000.0K)
    j_NH_Amax::typeof(50.0u"μmol*s^-1")
    j_NO_Amax::typeof(50.0u"μmol*s^-1")
    K_NH::typeof(10.0u"mmol*L^-1")
    K_NO::typeof(10.0u"mmol*L^-1")
    K_H::typeof(10.0u"mol*L^-1")
    X_NH::typeof(5.0u"mmol*L^-1")
    X_NO::typeof(10.0u"mmol*L^-1")
    X_H::typeof(10.0u"mol*L^-1")
    ρNO::typeof(0.7u"mmol*L^-1")
    y_E_CH_NH::typeof(1.25u"mol*mol^-1")
    y_EC_ECT::typeof(1.0u"mol*mol^-1")
end

function deb_settings(tspan; environment=[], use_environment=false, save_intermediate=false, timestep_days=1.0/24.0)
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

    load_settings(DEBSettings, DEBPars, OrderedDict(
        :settings => OrderedDict(
            :u0 => u0,
            :tspan => tspan,
            :environment => environment,
            :use_environment => use_environment,
            :apply_environment! => f() = nothing,
            :save_intermediate => save_intermediate,
            :timestep_days => timestep_days,
            :state_type => state_type
        ),
        :structures => OrderedDict(
            :Leaf => OrderedDict(
                :functions => OrderedDict(
                    :area => area_mass_kooijman,
                    :assim => shoot_assimilation!,
                    :assim_sub => photosynthesis_kooijman,
                    :rate => find_rate,
                ),
                :params => OrderedDict(
                    :J_L_F => (watts_to_light_mol(800.0)u"mol*m^2*s^-1", "mol/m²s, flux of useful photons"),
                    # TODO: should this acclimatise over time? this would require another state variable.

                    # J_L_F should feed back to affect plant state: photosynthesis as a sensory
                    # as well as energetic process [_@huner1998energy
                    :J_L_K => (watts_to_light_mol(300.0)u"mol*m^2*s^-1", "mol/m²s, half-saturation flux of useful photons"),
                    # Max specific uptake parameters that relate uptake to active surface area
                    # These are usually measured in μmol m⁻² s⁻¹ (wikipedia:Photosynthetic capacity).
                    # From [_@walker2014relationship :...full range of photosynthetically active radiation
                    # (PAR 0–1500 μmol·m−2·s−1) and three levels of Vcmax (25, 50 & 90 μmol·m−2·s−1)
                    :j_L_Amax => (20.0u"μmol*m^2*s^-1", "umol.m⁻².s⁻¹, max spec uptake of useful photons"),
                    :j_C_Amax => (90.0u"μmol*m^2*s^-1", "mol.m⁻².s⁻¹, max spec uptake of carbon dioxide"),
                    :j_O_Amax => (0.001u"μmol*m^2*s^-1", "mol.m⁻².s⁻¹, max spec uptake of oxygen"),
                    # Binding rates of gases to quantify photo-respiration
                    # Is this braodly similar accross plants?
                    :k_C_binding => (1.0u"mol*s^-1", :temp, "mols.s⁻¹, scaling rate for carbon dioxide"),
                    :k_O_binding => (1.0u"mol*s^-1", :temp, "mols.s⁻¹, scaling rate for oxygen"),
                    # Turnover. This controls metabolic rate along with area/mass relation and
                    # current reserves.
                    :k_E => (0.2u"mol*d^-1", (0.0u"mol*d^-1", 1.0u"mol*d^-1"), (:exposed, :time, :temp), "shoots reserve turnover rate"),
                    :k_EC => (0.2u"mol*d^-1", (0.0u"mol*d^-1", 1.0u"mol*d^-1"), (:exposed, :time, :temp), "shoots C-reserve turnover rate"),
                    :k_EN => (0.2u"mol*d^-1", (0.0u"mol*d^-1", 1.0u"mol*d^-1"), (:exposed, :time, :temp), "shoots N-reserve turnover rate"),
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
                    :κEN => (0.5, (0.0, 1.0), (:exposed, :time), "shoots non-processed N-reserve returned to N-reserve"),
                    :κsoma => (0.6, (0.0, 1.0), (:exposed,), "shoots reserve flux allocated to soma"),
                    :κrep => (0.05, (0.0, 1.0), (:exposed,), "shoots reserve flux allocated to development/reprod."),
                    :y_V_E => (0.7u"mol*mol^-1", (:exposed, :time), "from shoots reserve to structure: 0.3 lost as CO2?"),
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
                    :n_N_EC => (0.0u"mol*mol^-1", "N/C in C-reserve"),
                    :n_N_EN => (10.0u"mol*mol^-1", "N/C in N-reserve"),
                    :n_N_E => (0.2u"mol*mol^-1", "N/C in reserve. This should be calculated, not constant. 1.8181??? (10/11 * 0.3)/1.5"),

                    # Parameters that link moles to grams (wet weight)
                    # These should be calculated from first principles, not preset.
                    :w_P => (25.0u"g*mol^-1", "mol-weight of shoot product (wood): 1mol * 12g/mol => 12g"),
                    :w_V => (25.0u"g*mol^-1", "mol-weight of shoot structure: 0.15mol * 14g/mol + 1mol * 12g/mol => 14.1g"),
                    :w_EC => (25.0u"g*mol^-1", "mol-weight of shoot C-reserve: 1mol * 12g/mol => 12g"),
                    :w_EN => (25.0u"g*mol^-1", "mol-weight of shoot N-reserve: 10mol * 14g/mol + 1mol * 12g/mol => 152g"),
                    :w_E => (25.0u"g*mol^-1", "mol-weight of shoot reserve: 0.2mol * 14g/mol + 1mol * 12g/mol => 14.8g"),

                    :REFERENCE_TEMP => (310.0K, (273.0K, 325.0K), "temp for which rate pars are given"),
                    :ARRH_TEMP => (2000.0K, (200.0K, 4000.0K), "Arrhenius temp"),
                    :LOWER_BOUNDARY => (280.0K, (273.0K, 325.0K), "lower boundary tolerance range"),
                    :ARRH_LOWER => (20000.0K, (2000.0K, 40000.0K), "Arrhenius temp for lower boundary"),
                    :UPPER_BOUNDARY => (315.0K, (273.0K, 325.0K), "upper boundary tolerance range"),
                    :ARRH_UPPER => (70000.0K, (7000.0K, 140000.0K), "Arrhenius temp for upper boundary"),
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

                    :κEC => (0.5, (0.0, 1.0), (:exposed,), "roots  non-processed C-reserve returned to C-reserve"),
                    # the remaining fraction is translocated to the shoot
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
