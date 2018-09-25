using Base.Test

include("setup.jl")
include("tests/test_parameters.jl")
include("assimilation_inputs.jl")
include("solvers.jl")

# The change of the reserve energy E in time t can be written as
# d/dt * E = ˙pA − ˙pC, 
# where the assumption is that the mobilisation rate of reserve, ˙pC, is some function of the amount of
# reserve energy E and of structural volume V only. This function is fully determined by the


@testset "flux" begin
    sps = STATES_PER_STRUCTURE
    tspan = (0.0, 100.0)
    params = test_parameters(tspan)
    model = ModelWrapper(plant_model!, params, 0.0)
    t = 1.0 
    u = deepcopy(M0)
    du = deepcopy(M0)
    model(t, u, du)
    s = params.structures
    Jbase = params.Jbase[:,:,1]
    J = model.J
    J1 = model.J1
    loss = params.J1base[:,los, 1]

    du[1:5] + du[6:10]
    loss[1:5] + loss[6:10]
    du
    loss

    @testset "C loss" begin
        loss = params.J1base[:,los, 1]
        @test sum(du) + sum(loss) == 0.0
    end

    @testset "N loss" begin
        loss = params.J1base[:,los, 1]
        Nloss = loss[EN] * s[1].params[_n_N_EN] + loss[E] * s[1].params[_n_N_E]
        Nloss += loss[EN+sps] * s[2].params[_n_N_EN] + loss[E+sps] * s[2].params[_n_N_E]
        Ndu = du[EN] * s[1].params[_n_N_EN] + du[E] * s[1].params[_n_N_E]  + du[V] * s[1].params[_n_N_V]
        Ndu += du[EN+sps] * s[2].params[_n_N_EN] + du[E+sps] * s[2].params[_n_N_E]  + du[V+sps] * s[2].params[_n_N_V]
        @test Ndu + Nloss == 0.0
    end

    collect(zip(STATE_NAMES, du))

    @testset "growth" begin
        for (i, str) in enumerate(params.structures)
            CN = (J[i][i][V,gro] * (1/str.params[_y_V_E]) + J[i][i][E,gro])
            @test J[i][i][EC,gro] ≈ -CN/str.params[_y_E_CH_NO]
            @test J[i][i][EN,gro] ≈ -CN/str.params[_y_E_EN]
        end
    end

    @testset "translocation" begin
        κtra1 = calc_κtra(s[1].params)
        κtra2 = calc_κtra(s[2].params)

        Edrain = -κtra1 * (J1[1][EE,cat] + J1[1][ECN,cat])
        Egain = s[2].params[_y_E_ET] * κtra2 * J1[2][E,cat]
        Jtra1 = Edrain + Egain 
        lostE1 = κtra1 * (J1[1][E,cat]) * (1-s[1].params[_y_E_ET])

        Edrain = -κtra2 * (J1[2][EE,cat] + J1[2][ECN,cat])
        Egain = s[1].params[_y_E_ET] * κtra1 * J1[1][E,cat]
        Jtra2 = Edrain + Egain 
        lostE2 = κtra2 * (J1[2][E,cat]) * (1-s[2].params[_y_E_ET])

        @test Jtra1 + Jtra2 + lostE1 + lostE2 ≈ 0.0 atol=1e-20
    end

end

