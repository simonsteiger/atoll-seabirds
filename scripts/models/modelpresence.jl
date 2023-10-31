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
benchmark = true

# Load saved chains?
load = false

# Save the result?
save = true
!save && @warn "Samples will NOT be saved automatically."

# If not loading a chain, save results to path below
chainpath = "predictpresence_s_in_n.jls"

### MODEL SPECIFICATION ###

lu(x) = length(unique(x))

@model function modelpresence(
    r, s, n, PC, y,
    idx_sn, u_n, u_sn, Nv, Ng, Nb;
    Nr=lu(r), Ns=lu(s), Nn=lu(n), NPC=size(PC, 2), idx_sr=idx(s, r)
)

    # Priors for species × region
    μ_sxr ~ Normal()
    τ_sxr ~ Exponential(1)
    z_sxr ~ filldist(Normal(), Ns * Nr)
    α_sxr = μ_sxr .+ τ_sxr .* z_sxr

    # Priors for nesting types × PCs
    μ_pxn ~ filldist(Normal(), Nn)
    τ_pxn ~ filldist(Exponential(1), Nn)
    z_pxb ~ filldist(Normal(), Nb, NPC)
    z_pxg ~ filldist(Normal(), Ng, NPC)
    z_pxv ~ filldist(Normal(), Nv, NPC)
    z_pxn = reduce(vcat, [z_pxb, z_pxg, z_pxv])
    β_pxn = getindex.((μ_pxn,), u_n) .+ getindex.((τ_pxn,), u_n) .* getindex.((z_pxn,), u_sn)

    # Likelihood
    v = α_sxr[idx_sr] + sum(β_pxn[idx_sn, :] .* PC, dims=2)
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
    presence,
    num_species_within_nesting,
    unique_nesting,
    unique_species_within_nesting,
    count_species_by_nesting...,
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

if load
    chain = deserialize("chains/$chainpath")
end;

θ = generated_quantities(model, chain)

function predictpresence(α, β, n, s, r, X; idx_sr=idx(s, r))
    [rand.(BernoulliLogit.(α[i][idx_sr] .+ sum(β[i][n, :] .* X, dims=2))) for i in eachindex(α)]
end

α = [θ[i].α_sxr for i in eachindex(θ[:, 1])]
β = [θ[i].β_pxn for i in eachindex(θ[:, 1])]

nnu_long = [fill(num_nesting_unknown, length(num_region_unknown))...;]
nru_long = [fill.(num_region_unknown, length(num_nesting_unknown))...;]
nsu_long = [fill(num_species_unknown, length(num_region_unknown))...;]
PCu_long = reduce(vcat, [permutedims(hcat(fill(s, length(num_nesting_unknown))...)) for s in eachslice(PC_unknown, dims=1)])

preds = Matrix{Float64}(reduce(hcat, vec(predictpresence(α, β, nnu_long, nsu_long, nru_long, PCu_long));))

pct_preds = vec(mean(preds, dims=2))

ssu_long = [fill(Preprocess.str_species_unknown, length(num_region_unknown))...;]

sau_long = [fill.(Preprocess.str_atoll_unknown, length(num_nesting_unknown))...;]

df_preds = @chain begin
    DataFrame([sau_long, ssu_long, pct_preds], [:atoll, :species, :percent])
    unstack(_, :species, :percent)
end

CSV.write("../../data/newpreds.csv", df_preds)

# Find out what pushes some species onto the maledives 
# Not regional parameters
# Model parameters for these nesting types also not crazy
# So there are some outliers in the observed PC matrix? 
# PC histograms?
# Compare the observed and predicted data - are they very different?
