# Prior predictive checks

using Turing, StatsPlots

# Looks like TDist(3) is a good regularising prior
@model function sample_prior(y, Dist)
    intercept ~ Dist

    pc1 ~ Dist
    pc2 ~ Dist
    pc3 ~ Dist
    pc4 ~ Dist
    pc5 ~ Dist
    pc6 ~ Dist

    v = logistic(intercept + pc1 + pc2 + pc3 + pc4 + pc5 + pc6)
    y ~ Bernoulli(v)
end;

function logitpred(chn)
    res = logistic.(chn[:intercept] .+ chn[:pc1] .+ chn[:pc2] .+ chn[:pc3] .+ chn[:pc4] .+ chn[:pc5] .+ chn[:pc6])
    return res
end

pT = sample(sample_prior(missing, TDist(3)), Prior(), 100_000)
pN = sample(sample_prior(missing, Normal(0, 10)), Prior(), 100_000)

p = density()
density!([logitpred(x) for x in [pT, pN]], label=["T" "N"])
