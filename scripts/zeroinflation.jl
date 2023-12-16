# Probabilistic programming
using Turing, ReverseDiff, ParetoSmooth
# Benchmarking
using TuringBenchmarking, BenchmarkTools
# Model speed optimization
using LazyArrays
# Statistics
using StatsFuns, LinearAlgebra
# Working with tabular data
using Chain, DataFrames
# Plotting
using StatsPlots
# Saving results and logging
using Serialization, CSV, Dates, Markdown
# Random seeds
using Random
# New methods for these
import Distributions: logpdf, rand

# Path
const ROOT = dirname(Base.active_project())

# Load custom modules
include("$ROOT/src/count.jl")
include("$ROOT/src/utilities.jl")

# Make custom modules available
using .CountVariables
using .CustomUtilityFuns

struct ZIPoisson{T1,T2} <: DiscreteUnivariateDistribution
    λ::T1
    w::T2
end

function logpdf(d::ZIPoisson, y::Int)
    if y == 0
        logsumexp([log(d.w), log(1 - d.w) - d.λ])
    else
        log(1 - d.w) + logpdf(Poisson(d.λ), y)
    end
end

function rand(d::ZIPoisson)
    rand() <= d.w ? 0 : rand(Poisson(d.λ))
end

rand(d::ZIPoisson, N::Int) = map(_ -> rand(d), 1:N)

lun(x) = length(unique(x))

@model function zipmodel(
    r, s, n, PC,
    idx_sn, u_n, u_sn, Nv, Ng, Nb,
    y;
    Nr=lun(r), Ns=lun(s), Nn=lun(n), NPC=size(PC, 2), idx_sr=idx(s, r)
)

    # Priors for species × region
    μ_p_sxr ~ filldist(Normal(-1.5, 0.1), Ns)
    τ_p_sxr ~ filldist(InverseGamma(3, 0.5), Ns)
    z_p_sxr ~ filldist(Normal(), Ns, Nr)
    α_p_sxr = μ_p_sxr .+ τ_p_sxr .* z_p_sxr
    
    μ_λ_sxr ~ filldist(Normal(1, 0.1), Ns)
    τ_λ_sxr ~ filldist(InverseGamma(3, 0.1), Ns)
    z_λ_sxr ~ filldist(Normal(), Ns, Nr)
    α_λ_sxr = μ_λ_sxr .+ τ_λ_sxr .* z_λ_sxr

    # Priors for nesting types × PCs
    μ_p_pxn ~ filldist(Normal(0, 0.1), Nn, NPC)
    τ_p_pxn ~ filldist(InverseGamma(3, 0.1), Nn, NPC)
    z_p_pxb ~ filldist(Normal(), Nb, NPC)
    z_p_pxg ~ filldist(Normal(), Ng, NPC)
    z_p_pxv ~ filldist(Normal(), Nv, NPC)
    z_p_pxn = ApplyArray(vcat, z_p_pxb, z_p_pxg, z_p_pxv)
    β_p_pxn = μ_p_pxn[u_n, :] .+ τ_p_pxn[u_n, :] .* z_p_pxn[u_sn, :]
    
    μ_λ_pxn ~ filldist(Normal(1, 0.1), Nn, NPC)
    τ_λ_pxn ~ filldist(InverseGamma(3, 0.1), Nn, NPC)
    z_λ_pxb ~ filldist(Normal(), Nb, NPC)
    z_λ_pxg ~ filldist(Normal(), Ng, NPC)
    z_λ_pxv ~ filldist(Normal(), Nv, NPC)
    z_λ_pxn = ApplyArray(vcat, z_λ_pxb, z_λ_pxg, z_λ_pxv)
    β_λ_pxn = μ_λ_pxn[u_n, :] .+ τ_λ_pxn[u_n, :] .* z_λ_pxn[u_sn, :]

    # Likelihood
    p = logistic.(α_p_sxr[idx_sr] + sum(β_p_pxn[idx_sn, :] .* PC, dims=2))
    λ = exp.(α_λ_sxr[idx_sr] + sum(β_λ_pxn[idx_sn, :] .* PC, dims=2))
    y .~ ZIPoisson.(λ, p)

    # Generated quantities
    return (; y, p, λ, α_p_sxr, α_λ_sxr, β_p_pxn, β_λ_pxn)
end;

function simulate(params, r, s, n, X, idx_sn, u_n, u_sn, Nv, Ng, Nb; idx_sr=idx(s, r))
    map(params) do param
        p = logistic.(param.α_p_sxr[idx_sr] + sum(param.β_p_pxn[idx_sn, :] .* X, dims=2))
        λ = exp.(param.α_λ_sxr[idx_sr] + sum(param.β_λ_pxn[idx_sn, :] .* X, dims=2))
        rand.(ZIPoisson.(λ, p))
    end
end

inputs = [
    num_region_known,
    num_species_known,
    num_nesting_known,
    PC_known,
    num_species_within_nesting_known,
    unique_nesting_known,
    unique_species_within_nesting_known,
    count_species_by_nesting...,
]

counts = Int64.(nbirds)

model = zipmodel(inputs..., counts);

chain_prior = sample(model, Prior(), 1000);

θ_prior = peaceful_generated_quantities(model, chain_prior);

let k = keys(θ_prior[1])
    params = vec(getsamples(θ_prior, k...))
    priorsamples = reduce(vcat, simulate(params, inputs...))
    density(priorsamples, fillrange=0, fillalpha=0.2, legend=false)
end

backends = [
    Turing.Essential.ForwardDiffAD{0}(),
    Turing.Essential.ReverseDiffAD{false}(),
    Turing.Essential.ReverseDiffAD{true}()
];

TuringBenchmarking.run(TuringBenchmarking.make_turing_suite(model, adbackends=backends);)

# We cannot cache the model because of runtime if-statements (in rand() and logpdf() methods for ZIPoisson)

Turing.setadbackend(:reversediff)

# Configure sampling
sampler = NUTS(1000, 0.95; max_depth=10)
nsamples = 500
nthreads = 4
ndiscard = 50

chain = sample(model, sampler, nsamples; discard_initial=ndiscard)

serialize("$ROOT/results/chains/ziptest.jls", chain)