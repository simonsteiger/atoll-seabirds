using CSV, Serialization, DataFrames

# Pathing
const ROOT = dirname(Base.active_project())

include("$ROOT/scripts/visualization/diagnosticplots.jl")
using .DiagnosticPlots

chain_narrow = deserialize("$ROOT/scripts/models/chains/presence_narrow.jls")
chain_wide = deserialize("$ROOT/scripts/models/chains/presence_wide.jls")
chain_default = deserialize("$ROOT/scripts/models/chains/presence_default.jls")

preds_narrow = CSV.read("$ROOT/data/presencepreds_narrow.csv", DataFrame)
preds_wide = CSV.read("$ROOT/data/presencepreds_wide.csv", DataFrame)
preds_default = CSV.read("$ROOT/data/presencepreds_default.csv", DataFrame)

using StatsPlots

heatmap(Matrix(preds_default[!, 2:end] .- preds_narrow[!, 2:end]))

heatmap(Matrix(preds_default[!, 2:end] .- preds_wide[!, 2:end]))

using Chain

ps = map([0.7:0.05:0.9;]) do threshold
    @chain begin
        (preds_narrow[!, 2:end] .> threshold) .- (preds_wide[!, 2:end] .> threshold)
        Matrix{Float64}(_)
        heatmap(_, title="$threshold: $(Int64(sum(_))) species changed")
    end
end

plot(ps..., layout=(2, 3), size=(1200, 800))

map([0.7:0.05:0.95;]) do threshold
    x = sum(Matrix(preds_default[!, 2:end] .> threshold))
    "$threshold: $(Int64(x)) positive"
end

plot(
    map([0.7:0.05:0.9;]) do threshold
        heatmap(Matrix(preds_default[!, 2:end] .> threshold), title="$threshold")
    end...,
    layout=(2, 3),
    size=(1200, 800)
)