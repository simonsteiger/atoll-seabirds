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
chainpath = "predictpresence_new2.jls"
modelpath = "newmodel.jls"

chain = deserialize("chains/$chainpath")
model = deserialize("chains/$modelpath")

predictpresence(α, β, n, X) = rand.(Bernoulli.(α .+ β[n, :] * X))

θ = generated_quantities(model, chain)

α = [θ[i].α_sxr for i in eachindex(θ[:, 1])]
β = [θ[i].β_pxn for i in eachindex(θ[:, 1])]

long_X = [fill(PC_unknown, length(num_nesting_unknown))...;]

predictpresence.(α, β, num_nesting_unknown, )