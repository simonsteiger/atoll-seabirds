module CustomStatsFuns

export standardise,
       unstandardise

using StatsBase

# Z standardise
function standardise(x::AbstractArray)
    x̄, σ = mean(x), std(x)
    [(x - x̄) / σ for x in x]
end

# Recover natural scale from Z values
unstandardise(x::AbstractArray, x̄::Real, σ::Real) = x * σ .+ x̄

end