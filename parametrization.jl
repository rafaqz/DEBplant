using DiffEqBayes, ApproxBayes

tstop = 24 * 20
prob = DiscreteProblem(model, u, (1, tstop))
alg = FunctionMap(scale_by_time = true)
priors = metaflatten(Vector, model, FieldMetadata.prior)

abc_inference(prob::DEProblem, alg, t, data, priors; ϵ=0.001,
     distancefunction = euclidean, ABCalgorithm = ABCSMC, progress = false,
     num_samples = 500, maxiterations = 10^5, kwargs...)

