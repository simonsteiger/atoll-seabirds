# Full model to predict species presence

const PATH = "scripts/models"

# Paths relative to this folder
cd(PATH)

### PACKAGES AND MODULES ###

# Probabilistic programming
using Turing, TuringBenchmarking, LazyArrays
# Statistics
using StatsFuns, LinearAlgebra
# Working with tabular data
using Chain, DataFrames
# Plotting
using StatsPlots
# Saving results and logging
using Serialization, CSV, Dates

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

# Benchmark model?
benchmark = false

# Save the result?
save = true

# If not loading a chain, save results to path below
filename = "predictpresence_new.jls"

save && @info "The samples will be saved in $(joinpath(PATH, filename))."
!save && @info "The samples will NOT be saved automatically."

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
    α_sxr = μ_sxr .+ τ_sxr .* z_sxr[idx_sr]

    # Priors for nesting types × PCs
    μ_pxn ~ Normal()
    τ_pxn ~ Exponential(1)
    z_pxn ~ filldist(Normal(), Nn, NPC)
    β_pxn = μ_pxn .+ τ_pxn .* z_pxn

    # Likelihood
    y ~ arraydist(LazyArray(@~ BernoulliLogit.(α_sxr + β_pxn[n, :] * PC)))

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

@info "Starting sampling: $(Dates.now())"
chain = sample(model, sampler, MCMCThreads(), nsamples, nthreads; discard_initial=ndiscard);
@info "Finished sampling: $(Dates.now())"

save && serialize("chains/$filename", chain)

# Check acceptance rate
# plot_diagnostic(chain, "rhat")
# 
# θs = ["θ$i$j" for i in 1:6, j in 1:3]
# 
# params = extractparams(chain, ["λ", θs...])
# 
# preds = Dict()
# 
# for i in eachindex(num_species_unknown)
#     preds[Preprocess.str_species_unknown[i]] =
#         prediction(
#             "presence",
#             params,
#             num_species_unknown[i],
#             num_species_within_nesting_unknown[i],
#             num_nesting_unknown[i],
#             num_region_unknown,
#             num_atoll_unknown,
#             PC_unknown
#         )
# end
# 
# for k in keys(preds)
#     preds[k] = mean.(preds[k])
# end
# 
# df = DataFrame([preds...])
# insertcols!(df, 1, :atoll => Preprocess.str_atoll_unknown)
# 
# CSV.write("../../data/$savetofile.csv", df)
# 
# function specieshisto(df, idx)
#     p = histogram(df[:, idx], bins=10)
#     xlims!(0, 1)
#     title!(names(df)[idx])
#     return p
# end
# 
# specieshisto(df, 37)
# 
# TODO
# Add human population? PC or directly?
# Try TDist
# Prior predictive simulation - how sensitive is our model to variation in the predictors? Shouldn't be suuuper sensitive
# LOO - how sensitive is our model to single data points in the known data
