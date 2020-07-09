using PrettyTables, Latexify

dir = @__DIR__
include(joinpath(dir, "load.jl"))

fn = fieldnameflatten(model)
val = flatten(model)
remdups(a) = begin
    lastc = ""
    for (i, c) in enumerate(a)
        if c == lastc
            a[i] = ""
        else
            lastc = c
        end
    end
    a
end
comp = remdups([string.(parenttypeflatten(model))...])

desc = metaflatten(model, FieldMetadata.description)

data = hcat([comp...], [fn...], [val...], [desc...])
headers = ["Component", "Parameter", "Value", "Description"]

pretty_table(data, headers, markdown)



q = :((J_C_M, J_N_M) = j_E_mai * tempcorr * V / (y_E_EC,  y_E_EN))

latexify(q)
