using Revise

using Base.Test
using RCall
using BiophysicalModels
using DynamicEnergyBudgets
using DynamicEnergyBudgets.Flux
using DynamicEnergyBudgets.build_axis
using DynamicEnergyBudgets.split_flux
using DynamicEnergyBudgets.TRANS
using DataStructures

include("/home/raf/Uni/Masters/julia/DynamicEnergyBudgets/test/test_settings.jl")

const BI_XTOL = 1e-15

state_names = fieldnames(StatePVCNE)
u0 = [0.0, 1e-4, 1e-4, 1e-4, 1e-4, 0.0, 1e-4, 1e-4, 1e-4, 10.0]

function transform_R_to_jul(Jr)
    Jr = reshape(Jr, (10, 8))
    Jbase1 = build_axis(state_names, TRANS, 1, 1:1)
    Jbase2 = build_axis(state_names, TRANS, 1, 1:1)
    Jjul = (split_flux(Jbase1, 1), split_flux(Jbase2, 1))
    Jjul[1][:, :ass] = Jr[1:5, 1]
    Jjul[1][:, :gro] = Jr[1:5, 2]
    Jjul[1][:, :mai] = Jr[1:5, 3]

    Jjul[2][:, :ass] = Jr[1:5, 1 + 5]
    Jjul[2][:, :gro] = Jr[6:10, 2 + 5]
    Jjul[2][:, :mai] = Jr[6:10, 3 + 5]

    Jjul[1][1:5, :rep] = Jr[1:5, 4]
    Jjul[1][1:5, :tra] = Jr[1:5, 5]
    Jjul[2][1:5, :rep] = Jr[6:10, 4]
    Jjul[2][1:5, :tra] = Jr[6:10, 5]
    return Jjul
end

function test_states(name, M_ref, M_test, tol)
    for i = 1:10
        r = M_ref[i]
        t = M_test[i]
        println(string(name, ": M_", state_names[rem(i-1,5) + 1], " r = t : ", r, " = ", t))
        @test r ≈ t atol = 1e-15
    end
end

function test_fluxes(name, J_ref, J_test)
    J_test[2][:C,:mai] += J_test[2][:C,:rej]
    J_test[1][:N,:mai] += J_test[1][:N,:rej]
    for i = 1:2
        for k = 1:5
            for l in eachindex(TRANS)
                if TRANS[l] == :rej break end 
                r = J_ref[i][k, l]
                t = J_test[i][k, l]
                println(string(name, ": J[", i, "][", state_names[k], ",", 
                               DynamicEnergyBudgets.TRANS[l], "] r = t : ", r, " = ", t))
                @test r ≈ t atol = 1e-15
            end
        end
    end
end

@testset "Metabolic rate is correct" begin
    tspan = (0.0, 10.0)
    t = 0.0
    u = deepcopy(u0)
    du = deepcopy(u0)
    settings = test_settings(tspan)
    BiophysicalModels.runmodel!(du, settings, u, t)
    rS = settings.structures[1].rates[1]
    rR = settings.structures[2].rates[1]
    @test rS ≈ 0.0756394041835524
    @test rR ≈ 0.1999960679199535296124
end

@testset "runmodel! integrates a single timestep" begin
    tspan = (0.0, 100.0)
    settings = test_settings(tspan)
    u = deepcopy(u0)
    du = deepcopy(u0)
    t = 1.0 
    BiophysicalModels.runmodel!(du, settings, u, t)

    @rput u0
    R"""
    library(deSolve)
    setwd('/home/raf/Uni/Masters/julia/DynamicEnergyBudgets/test/original')
    source('findr.R')
    source('flux_plant.R')
    source('pars_plant.R')
    source('tempcorr.R')
    options(digits=22)
    Out1 <- flux_plant(0, u0, p)
    JM1 = Out1$J1
    J1 <- Out1$J2
    setwd('..')
    """
    @rget J1 JM1

    J1_ref = transform_R_to_jul(J1)
    J1_test = (settings.structures[1].J, settings.structures[2].J)
    test_fluxes("J1", J1_ref, J1_test)
end

function get_J(structures, n::Int)
    (structures[1].Jbase[:,:,n], structures[2].Jbase[:,:,n])
end

@testset "Short timespan: no assimilation" begin

    tspan = (0.0, 10.0)
    settings = test_settings(tspan; save_intermediate=true)
    s = settings.structures
    sol = integrate(settings)
    M9 = sol.u[10]

    @rput u0 M9
    R"""
    library(deSolve)
    setwd('/home/raf/Uni/Masters/julia/DynamicEnergyBudgets/test/original')
    source('findr.R')
    source('flux_plant.R')
    source('pars_plant.R')
    source('tempcorr.R')
    options(digits=22)
    Out1 <- flux_plant(0, u0, p)
    JM1 = Out1$J1
    J1 <- Out1$J2
    Out10 <- flux_plant(0, M9, p)
    JM10 = Out10$J1
    J10 <- Out10$J2
    setwd('..')
    """
    @rget J1 JM1 J10 JM10

    J1_ref = transform_R_to_jul(J1)
    J10_ref = transform_R_to_jul(J10)
    J1_test = get_J(s, 2)
    J10_test = get_J(s, 11)
    test_fluxes("J1", J1_ref, J1_test)
    test_fluxes("J10", J10_ref, J10_test)

    M1_ref = u0 + JM1
    M1_test = sol.u[2]
    M10_ref = M9 + JM10
    M10_test = sol.u[11]
    test_states("M1", M1_ref, M1_test, 1e-15)
    test_states("M10", M10_ref, M10_test, 1e-15)
end

@testset "Medium timespan: assimilation active" begin
    tspan = (0.0, 60.0)
    settings = test_settings(tspan; save_intermediate=true)
    s = settings.structures
    sol = integrate(settings)
    M19 = sol.u[20]
    M39 = sol.u[40]
    M59 = sol.u[60]

    @rput u0 M19 M39 M59
    R"""
    library(deSolve)
    setwd('/home/raf/Uni/Masters/julia/DynamicEnergyBudgets/test/original')
    source('findr.R')
    source('flux_plant.R')
    source('pars_plant.R')
    source('tempcorr.R')
    options(digits=22)
    Out20 <- flux_plant(0, M19, p)
    JM20 = Out20$J1
    J20 <- Out20$J2
    Out40 <- flux_plant(0, M39, p)
    JM40 = Out40$J1
    J40 <- Out40$J2
    Out60 <- flux_plant(0, M59, p)
    JM60 = Out60$J1
    J60 <- Out60$J2
    setwd('..')
    """
    @rget J1 JM1 J10 JM10 J20 JM20 J40 JM40 J60 JM60

    J40_ref = transform_R_to_jul(J40)
    J60_ref = transform_R_to_jul(J60)
    J40_test = get_J(s, 41)
    J60_test = get_J(s, 61)
    test_fluxes("J40", J40_ref, J40_test)
    test_fluxes("J60", J60_ref, J60_test)

    M20_ref = M19 + JM20
    M40_ref = M39 + JM40
    M60_ref = M59 + JM60
    M20_test = sol.u[21]
    M40_test = sol.u[41]
    M60_test = sol.u[61]

    test_states("M20", M20_ref, M20_test, 1e-15)
    test_states("M40", M40_ref, M40_test, 1e-15)
    test_states("M60", M60_ref, M60_test, 1e-15)
end

@testset "Large timespan: should have everything happening" begin
    tspan = (0.0, 500.0)
    settings = test_settings(tspan; save_intermediate=true)
    s = settings.structures
    sol = integrate(settings)
    M69 = sol.u[70]
    M79 = sol.u[80]
    M99 = sol.u[100]
    M199 = sol.u[200]
    M499 = sol.u[500]

    @rput M69 M79 M99 M199 M499
    R"""
    library(deSolve)
    setwd('/home/raf/Uni/Masters/julia/DynamicEnergyBudgets/test/original')
    source('findr.R')
    source('flux_plant.R')
    source('pars_plant.R')
    source('tempcorr.R')
    options(digits=22)
    Out70 <- flux_plant(0, M69, p)
    JM70 = Out70$J1
    J70 <- Out70$J2
    Out80 <- flux_plant(0, M79, p)
    JM80 = Out80$J1
    J80 <- Out80$J2
    Out100 <- flux_plant(0, M99, p)
    JM100 = Out100$J1
    J100 <- Out100$J2
    Out200 <- flux_plant(0, M199, p)
    JM200 = Out200$J1
    J200 <- Out200$J2
    Out500 <- flux_plant(0, M499, p)
    JM500 = Out500$J1
    J500 <- Out500$J2
    setwd('..')
    """
    @rget J70 JM70 J80 JM80 J100 JM100 J200 JM200 J500 JM500

    J80_ref = transform_R_to_jul(J80)
    J100_ref = transform_R_to_jul(J100)
    J500_ref = transform_R_to_jul(J500)

    J80_test = get_J(s, 81)
    J100_test = get_J(s, 101)
    J500_test = get_J(s, 501)

    test_fluxes("J80", J80_ref, J80_test)
    test_fluxes("J100", J100_ref, J100_test)
    test_fluxes("J500", J500_ref, J500_test)

    M70_ref = M69 + JM70
    M80_ref = M79 + JM80
    M100_ref = M99 + JM100
    M200_ref = M199  + JM200
    M500_ref = M499  + JM500

    M70_test = sol.u[71]
    M80_test = sol.u[81]
    M100_test = sol.u[101]
    M200_test = sol.u[201]
    M500_test = sol.u[501]

    test_states("M70", M70_ref, M70_test, 1e-15)
    test_states("M80", M80_ref, M80_test, 1e-15)
    test_states("M100", M100_ref, M100_test, 1e-15)
    test_states("M200", M200_ref, M200_test, 1e-15)
    test_states("M500", M500_ref, M500_test, 1e-15)
end
