using JLD2, Microclimate, Unitful, DataStructures

loadenvironments(dir) = begin
    locationspath = joinpath(dir, "microclimate/locations.jld")
    @load locationspath tas desert qld
    environments = OrderedDict(:tas=>tas, :desert=>desert, :qld=>qld)
    env = first(values(environments))
    tspan = (0:1:length(radiation(env)) - 1) * hr
    environments, tspan
end
