using NicheMap

long = 147.29
lat = -36.87
nm = nichemap_global([long,lat]; years=10)
mountain = extract_microclim(nm)
