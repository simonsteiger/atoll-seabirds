module CustomStatsFuns

export standardise,
       unstandardise

using StatsBase

# Z standardise
standardise(x) = (x .- mean(x)) ./ std(x)

# Recover natural scale from Z values
unstandardise(x::AbstractArray, x̄::Real, σ::Real) = x * σ .+ x̄

end