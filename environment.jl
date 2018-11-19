function load_environment()
    # using IndexedTables
    environment = load(joinpath(dir, "environment.jld"))["environment"]
    # environment = nichemap_global("Adelaide", years=10)
    Microclimate.MicroclimateTable(
      Table(environment.soil),
      Table(environment.shadsoil),
      Table(environment.metout),
      Table(environment.shadmet),
      Table(environment.soilmoist),
      Table(environment.shadmoist),
      Table(environment.humid),
      Table(environment.shadhumid),
      Table(environment.soilpot),
      Table(environment.shadpot),
      Table(environment.plant),
      Table(environment.shadplant),
      environment.RAINFALL,
      environment.dim,
      environment.ALTT,
      environment.REFL,
      environment.MAXSHADES,
      environment.longlat,
      environment.nyears,
      environment.timeinterval,
      environment.minshade,
      environment.maxshade,
      environment.DEP,
    )
end
