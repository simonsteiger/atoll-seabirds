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

@info "packages loaded"

# Load custom modules
include("../preprocessing/countvars.jl")
include("../../src/postprocess.jl")
include("../../src/utilities.jl")
include("../visualization/diagnosticplots.jl")
include("../visualization/paramplots.jl")

@info "modules loaded"

# Make custom modules available
using .CountVariables
using .Postprocess
using .CustomUtilityFuns
using .DiagnosticPlots
using .ParamPlots

# Set seed
Random.seed!(42)

# Benchmark model?
benchmark = false

# Load saved chains?
load = true

# Save the result?
save = true
!save && @warn "Samples will NOT be saved automatically."

# If not loading a chain, save results to path below
chainpath = "chains_count.jls"

### MODEL SPECIFICATION ###

lu(x) = length(unique(x))

@model function modelcount(
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
    μ_pxn ~ filldist(Normal(), Nn, NPC)
    τ_pxn ~ filldist(Exponential(1), Nn, NPC)
    z_pxb ~ filldist(Normal(), Nb, NPC)
    z_pxg ~ filldist(Normal(), Ng, NPC)
    z_pxv ~ filldist(Normal(), Nv, NPC)
    z_pxn = ApplyArray(vcat, z_pxb, z_pxg, z_pxv)
    β_pxn = μ_pxn[u_n, :] .+ τ_pxn[u_n, :] .* z_pxn[u_sn, :]

    σ ~ Exponential(1)

    # Likelihood
    μ = vec(α_sxr[idx_sr] + sum(β_pxn[idx_sn, :] .* PC, dims=2))
    y ~ MvNormal(μ, σ^2 * I) # Can we LazyArray μ?

    # Generated quantities
    return (α_sxr=α_sxr, β_pxn=β_pxn)
end;

# Create model
model = modelcount(
    num_region_known,
    num_species_known,
    num_nesting_known,
    PC_known,
    log.(nbirds),
    num_species_within_nesting_known,
    unique_nesting_known,
    unique_species_within_nesting_known,
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
    sampler = NUTS(1000, 0.95; max_depth=10)
    nsamples = 20_000
    nthreads = 4
    ndiscard = 5000

    @info """Sampler: $(string(sampler))
    Samples: $(nsamples)
    Threads: $(nthreads)
    Discard: $(ndiscard)
    """

    @info "🚀 Starting sampling: $(Dates.now())"
    chain = sample(model, sampler, MCMCThreads(), nsamples, nthreads; discard_initial=ndiscard)

    save && serialize("chains/$chainpath", chain)
    isfile("chains/$chainpath") && @info "💾 Chain saved to '$(PATH)chains/$chainpath'."
end;

@info "chain loaded"

θ = generated_quantities(model, chain);

@info "parameters extracted"

function predictcount(α, β, σ, idx_sn, s, r, X; idx_sr=idx(s, r))
    out = Vector(undef, length(α))
    for i in eachindex(α)
        out[i] = rand.(Normal.(α[i][idx_sr] .+ sum(β[i][idx_sn, :] .* X, dims=2), σ[i]))
    end
    return out
end

function getsamples(θ, sym)
    @chain begin
        map(1:size(θ, 2)) do j
            [θ[i][sym] for i in eachindex(θ[:, j])]
        end
        reduce(vcat, _)
    end
end

α, β = [getsamples(θ, s) for s in [:α_sxr, :β_pxn]];
σ = reduce(hcat, get_params(chain).σ)'

threshold = ppres .> 0.8

countpreds_unknown = @chain begin
    predictcount(
        α,
        β,
        σ,
        num_species_within_nesting_unknown[threshold],
        num_species_unknown[threshold],
        num_region_unknown[threshold],
        PC_unknown[threshold, :],
    )
    reduce(hcat, _)
    Matrix{Float64}(_)
end

countpreds_known = @chain begin
    predictcount(
        α,
        β,
        σ,
        num_species_within_nesting_known,
        num_species_known,
        num_region_known,
        PC_known,
    )
    reduce(hcat, _)
    Matrix{Float64}(_)
end

avg_preds_unknown = vec(mean(countpreds_unknown, dims=2))
avg_preds_known = vec(mean(countpreds_known, dims=2))

pred_x = avg_preds_known[num_species_known .== 29]
obs_x = log.(nbirds)[num_species_known .== 29]
scatter(eachindex(pred_x), pred_x, label="pred")
scatter!(eachindex(pred_x), obs_x, label="obs")

sum(pred_x)
sum(obs_x)
# 
# ssu_long = [fill(Preprocess.str_species_unknown, length(num_region_unknown))...;]
# 
# sau_long = [fill.(Preprocess.str_atoll_unknown, length(num_nesting_unknown))...;]
# 
# df_preds = @chain begin
#     DataFrame([sau_long, ssu_long, pct_preds], [:atoll, :species, :percent])
#     unstack(_, :species, :percent)
# end
# 
# CSV.write("../../data/newpreds.csv", df_preds)
