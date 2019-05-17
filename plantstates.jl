u = zeros(length(STATE) * 2)mol
ulabelled = LArray{STATEKEYS}(u)

largeseed = copy(ulabelled)
largeseed.VS = 0.20mg / (25.0g/mol)
largeseed.CS = 5.0mg  / (25.0g/mol)
largeseed.NS = 0.20mg  / (25.0g/mol)
largeseed.VR = 0.05mg / (25.0g/mol)
largeseed.CR = 1.0mg  / (25.0g/mol)
largeseed.NR = 0.05mg  / (25.0g/mol)

smallseed = largeseed ./ 100

plant = copy(ulabelled)
plant.VS = 12.0g  / (25.0g/mol)
plant.CS = 4.0g   / (25.0g/mol)
plant.NS = 0.15g  / (25.0g/mol)
plant.VR = 4.0g   / (25.0g/mol)
plant.CR = 2.0g   / (25.0g/mol)
plant.NR = 0.1g  / (25.0g/mol)

# states = OrderedDict("Plant" => plant, "Large seed" => largeseed, "Small seed" => smallseed);
states = OrderedDict("Large seed" => largeseed);
nothing
