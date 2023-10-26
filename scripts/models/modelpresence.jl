# Full model to predict species presence

const PATH = "scripts/models/"

# Paths relative to this folder
cd(PATH)

### PACKAGES AND MODULES ###

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

# Add custom modules
using .Preprocess
using .Postprocess
using .CustomUtilityFuns
using .DiagnosticPlots
using .ParamPlots

### SETTINGS ###

Random.seed!(42)

# Benchmark model?
benchmark = true

# Save the result?
save = true
!save && @info "The samples will NOT be saved automatically."

# If not loading a chain, save results to path below
chainpath = "predictpresence_new2.jls"
modelpath = "newmodel.jls"

### MODEL SPECIFICATION ###

lu(x) = length(unique(x))
idx(i, j) = i .+ (j .- 1) * maximum(i)

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
    y ~ arraydist(LazyArray(@~ BernoulliLogit.(α_sxr[idx_sr] + β_pxn[n, :] * PC)))

    return (α_sxr=α_sxr, β_pxn=β_pxn)
end;

model = modelpresence(
    num_region,
    num_species,
    num_nesting,
    PC,
    presence
);

### BENCHMARK ###

# Benchmark different backends to find out which is fastest
let adbackends = [:forwarddiff, :reversediff, :reversediff_compiled]
    benchmark && benchmark_model(model; check=true, adbackends=adbackends)
end

### SAMPLING CONFIGURATION ###

# Set AD backend to :reversediff and compile with setrdcache(true)
Turing.setadbackend(:reversediff)
Turing.setrdcache(true)

sampler = NUTS(1000, 0.99; max_depth=10)
nsamples = 2000
nthreads = 3
ndiscard = 1000

@info """\nSampler: $(string(sampler))
Samples: $(nsamples)
Threads: $(nthreads)
Discard: $(ndiscard)
"""

### SAMPLE ###

@info "🚀 Starting sampling: $(Dates.now())"
chain = sample(model, sampler, MCMCThreads(), nsamples, nthreads; discard_initial=ndiscard);
@info "🏁 Finished sampling: $(Dates.now())"

if save
    serialize("chains/$chainpath", chain)
    serialize("chains/$modelpath", model)
end

if isfile("chains/$chainpath") && isfile("chains/$modelpath")
    @info """
    \n💾 Save successful!
    Chain saved to '$PATH$chainpath'.
    Model saved to '$PATH$modelpath'.
    """
end
