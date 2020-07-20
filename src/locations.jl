# Build location microclimates from Microclim datasets
# Mostly we just load the saved locations from JLD files

datapath = "/home/raf/Data/microclim_oz"
shade = 0
years = 2005:2011 # One extra year for overhang from 2010

function load_locations(datapath::String, years, shade)
    envgrid = load_grid(datapath, years, shade);

    t1 = MicroclimPoint(envgrid, CartesianIndex(65, 35))
    t2 = MicroclimPoint(envgrid, CartesianIndex(60, 35))
    r3 = MicroclimPoint(envgrid, CartesianIndex(55, 35))

    envgrid = nothing

    save(joinpath(dir, "data/locations.jld2"), Dict("t1" => t1, "t2" => t2, "t3" => t3))
end
