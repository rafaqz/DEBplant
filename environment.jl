using Revise
using TypedTables
using Microclimate
using DataFrames
using JLD2

function load_environment()
    @load "environment.jld"
    env
end

# using NicheMap
# env = nichemap_global([144, -37]; years=10) 
# env = convert(MicroclimateTables, env)
# @save "environment.jld" env
