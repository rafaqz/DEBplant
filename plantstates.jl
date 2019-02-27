u = zeros(12)mol
ulabelled = LArray{STATEKEYS}(u)

smallseed = copy(ulabelled)
smallseed.VS = 1e-3mg / (25.0g/mol)
smallseed.CS = 1e-3mg / (25.0g/mol)
smallseed.NS = 1e-4mg / (25.0g/mol)
smallseed.VR = 1e-3mg / (25.0g/mol)
smallseed.CR = 1.0mg  / (25.0g/mol)
smallseed.NR = 0.05mg / (25.0g/mol)

largeseed = copy(ulabelled)
largeseed.VS = 1e-1mg  / (25.0g/mol)
largeseed.CS = 1e-1mg  / (25.0g/mol)
largeseed.NS = 1e-1mg  / (25.0g/mol)
largeseed.VR = 1e-1mg  / (25.0g/mol)
largeseed.CR = 100.0mg / (25.0g/mol)
largeseed.NR = 5.0mg   / (25.0g/mol)

plant = copy(ulabelled)
plant.VS = 10g    / (25.0g/mol)
plant.CS = 10g    / (25.0g/mol)
plant.NS = 1.0g   / (25.0g/mol)
plant.VR = 5.0g   / (25.0g/mol)
plant.CR = 5.0g   / (25.0g/mol)
plant.NR = 0.5g   / (25.0g/mol)

states = Dict(:smallseed => smallseed, :largeseed => largeseed, :plant => plant);
nothing
