const PATH = "scripts/models/"

# Paths relative to this folder
cd(PATH)

# Probabilistic programming
using Turing, TuringBenchmarking, LazyArrays, Random
# Statistics
using StatsFuns, LinearAlgebra
# Working with tabular data
using Chain, DataFrames
# Plotting
using StatsPlots
# Saving results and logging
using Serialization, CSV, Dates, Markdown

# Load custom modules
include("../preprocessing/preprocess.jl")
include("../../src/postprocess.jl")
include("../../src/utilities.jl")
include("../visualization/diagnosticplots.jl")
include("../visualization/paramplots.jl")

# Make custom modules available
using .Preprocess
using .Postprocess
using .CustomUtilityFuns
using .DiagnosticPlots
using .ParamPlots

# Set seed
Random.seed!(42)

# Benchmark model?
benchmark = false

# Load saved chains?
load = false

# Save the result?
save = true
!save && @warn "Samples will NOT be saved automatically."

# If not loading a chain, save results to path below
chainpath = "predictpresence_new2.jls"

### MODEL SPECIFICATION ###

lu(x) = length(unique(x))

@model function modelpresence(r, s, n, PC, y;       # Main inputs
    Nr=lu(r), Ns=lu(s), Nn=lu(n), NPC=size(PC, 2),  # Number of groups
    idx_sr=idx(s, r))                               # Vector for indexing

    # Priors for species × region
    μ_sxr ~ Normal()
    τ_sxr ~ Exponential(1)
    z_sxr ~ filldist(Normal(), Ns * Nr)
    α_sxr = μ_sxr .+ τ_sxr .* z_sxr

    # Priors for nesting types × PCs
    μ_pxn ~ Normal()
    τ_pxn ~ Exponential(1)
    z_pxn ~ filldist(Normal(), Nn, NPC)
    β_pxn = μ_pxn .+ τ_pxn .* z_pxn

    # Likelihood
    v = α_sxr[idx_sr] + sum(β_pxn[n, :] .* PC, dims=2)
    y .~ BernoulliLogit.(v)

    # Generated quantities
    return (α_sxr=α_sxr, β_pxn=β_pxn)
end;

# Create model
model = modelpresence(
    num_region,
    num_species,
    num_nesting,
    PC,
    presence
);

# Benchmark different backends to find out which is fastest
let adbackends = [:forwarddiff, :reversediff, :reversediff_compiled]
    benchmark && benchmark_model(model; check=true, adbackends=adbackends)
end

# Sample from model unless a saved chain should be used
if load
    chain = deserialize("chains/$chainpath")
else
    # Set AD backend to :reversediff and compile with setrdcache(true)
    Turing.setadbackend(:reversediff)
    Turing.setrdcache(true)

    # Configure sampling
    sampler = NUTS(1000, 0.99; max_depth=10) # Sampler picks depth=9, restrict for speedup?
    nsamples = 2000
    nthreads = 3
    ndiscard = 1000

    @info """Sampler: $(string(sampler))
    Samples: $(nsamples)
    Threads: $(nthreads)
    Discard: $(ndiscard)
    """

    @info "🚀 Starting sampling: $(Dates.now())"
    chain = sample(model, sampler, MCMCThreads(), nsamples, nthreads; discard_initial=ndiscard)

    save && serialize("chains/$chainpath", chain)
    isfile("chains/$chainpath") && @info "💾 Chain saved to '$(PATH)chains/$chainpath'."
end

load && chain = deserialize("chains/$chainpath")

θ = generated_quantities(model, chain)

function predictpresence(α, β, n, s, r, X; idx_sr=idx(s, r))
    [rand.(BernoulliLogit.(α[i][idx_sr] .+ β[i][n, :] .* X)) for i in eachindex(α)]
end

α = [θ[i].α_sxr for i in eachindex(θ[:, 1])]
β = [θ[i].β_pxn for i in eachindex(θ[:, 1])]

nnu_long = [fill(num_nesting_unknown, length(num_region_unknown))...;]
nru_long = [fill(num_region_unknown, length(num_nesting_unknown))...;]
nsu_long = [fill(num_species_unknown, length(num_region_unknown))...;]
PCu_long = [fill(PC_unknown, length(num_nesting_unknown))...;]

preds = predictpresence(α, β, nnu_long, nsu_long, nru_long, PCu_long);
