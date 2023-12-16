module PresenceModels
    
using Turing, LazyArrays, StatsFuns
include("../src/utilities.jl")
using .CustomUtilityFuns

export mSp, mMp, mLp, mXLp, simulate

# Helper for model declaration
lu(x) = length(unique(x))

@model function mSp(r, s, y; Nr=lu(r), Ns=lu(s), idx_sr=idx(s, r))
    # Priors for species × region
    α_sxr ~ filldist(Beta(1, 2), Ns * Nr)

    # Likelihood
    y ~ arraydist(LazyArray(@~ Bernoulli.(α_sxr[idx_sr])))

    return (; α_sxr)
end;

function simulate(α, s, r; idx_sr=idx(s, r))
    [rand.(Bernoulli.(α[i][idx_sr])) for i in eachindex(α)]
end

@model function m2p(r, s, PC, y; Nr=lu(r), Ns=lu(s), NPC=size(PC, 2), idx_sr=idx(s, r))
    # Priors for species × region
    α_sxr ~ filldist(Beta(1, 2), Ns * Nr)

    # Priors for PCs, independent of nesting type
    β_p = filldist(Beta(1, 2), NPC)

    # Likelihood
    v = α_sxr[idx_sr] + sum(β_p .* PC, dims=2)
    y .~ Bernoulli.(v)
end;

@model function mMp(
    r, s, n, PC, # inputs
    idx_sn, u_n, u_sn, Nv, Ng, Nb, # indexes, lengths, etc
    y; # outcome
    Nr=lu(r), Ns=lu(s), Nn=lu(n), NPC=size(PC, 2), idx_sr=idx(s, r)
)

    # Priors for species × region
    α_sxr ~ filldist(Normal(0, 1), Ns * Nr)

    # Priors for nesting types × PCs
    μ_pxn ~ filldist(Normal(0, 0.2), Nn, NPC)
    σ_pxn ~ filldist(InverseGamma(3, 0.5), Nn, NPC)
    z_pxb ~ filldist(Normal(), Nb, NPC)
    z_pxg ~ filldist(Normal(), Ng, NPC)
    z_pxv ~ filldist(Normal(), Nv, NPC)
    z_pxn = ApplyArray(vcat, z_pxb, z_pxg, z_pxv)
    β_pxn = @. μ_pxn[u_n, :] + σ_pxn[u_n, :] * z_pxn[u_sn, :]

    # Likelihood
    v = α_sxr[idx_sr] + sum(β_pxn[idx_sn, :] .* PC, dims=2)
    y .~ BernoulliLogit.(v)

    # Generated quantities
    return (; α_sxr, β_pxn)
end;

# Same inputs as main model allows us to pass the same vector
function simulate(params, r, s, n, X, idx_sn, u_n, u_sn, Nv, Ng, Nb; idx_sr=idx(s, r))
    map(params) do p
        β = sum(p.β_pxn[idx_sn, :] .* X, dims=2)
        @. logistic(p.α_sxr[idx_sr] + β)
    end
end

@model function mLp(
    r, s, n, PC, # inputs
    idx_sn, u_n, u_sn, Nv, Ng, Nb, # indexes, lengths, etc
    y; # outcome
    Nr=lu(r), Ns=lu(s), Nn=lu(n), NPC=size(PC, 2), idx_sr=idx(s, r)
)

    # Species prior
    α_s ~ filldist(Normal(0, 0.5), Ns)
    # Region prior
    α_r ~ filldist(Normal(0, 0.5), Nr)
    # Priors for species × region
    α_sxr ~ filldist(Normal(0, 0.5), Ns * Nr)

    # Priors for nesting types × PCs
    μ_pxn ~ filldist(Normal(0, 0.2), Nn, NPC)
    σ_pxn ~ filldist(InverseGamma(3, 0.5), Nn, NPC)
    z_pxb ~ filldist(Normal(), Nb, NPC)
    z_pxg ~ filldist(Normal(), Ng, NPC)
    z_pxv ~ filldist(Normal(), Nv, NPC)
    z_pxn = ApplyArray(vcat, z_pxb, z_pxg, z_pxv)
    β_pxn = @. μ_pxn[u_n, :] + σ_pxn[u_n, :] * z_pxn[u_sn, :]

    # Likelihood
    β = sum(β_pxn[idx_sn, :] .* PC, dims=2)
    y ~ arraydist(LazyArray(@~ BernoulliLogit.(α_s[s] + α_r[r] + α_sxr[idx_sr] + β)))

    # Generated quantities
    return (; y, α_s, α_r, α_sxr, β_pxn)
end;

# Same inputs as main model allows us to pass the same vector
function simulate(params, r, s, n, X, idx_sn, u_n, u_sn, Nv, Ng, Nb; idx_sr=idx(s, r))
    map(params) do p
        β = sum(p.β_pxn[idx_sn, :] .* X, dims=2)
        @. logistic(p.α_s[s] + p.α_r[r] + p.α_sxr[idx_sr] + β)
    end
end

@model function mXLp(
    a, r, s, n, PC, # inputs
    idx_sn, u_n, u_sn, Nv, Ng, Nb, # indexes, lengths, etc
    y; # outcome
    Na=lu(a), Nr=lu(r), Ns=lu(s), Nn=lu(n), NPC=size(PC, 2), idx_sr=idx(s, r)
)

    # Atoll prior (rats?)
    α_a ~ filldist(Normal(0, 0.5), Na)
    # Species prior
    α_s ~ filldist(Normal(0, 0.5), Ns)
    # Region prior
    α_r ~ filldist(Normal(0, 0.5), Nr)
    # Priors for species × region
    α_sxr ~ filldist(Normal(0, 0.5), Ns * Nr)

    # Priors for nesting types × PCs
    μ_pxn ~ filldist(Normal(0, 0.2), Nn, NPC)
    σ_pxn ~ filldist(InverseGamma(3, 0.5), Nn, NPC)
    z_pxb ~ filldist(Normal(), Nb, NPC)
    z_pxg ~ filldist(Normal(), Ng, NPC)
    z_pxv ~ filldist(Normal(), Nv, NPC)
    z_pxn = ApplyArray(vcat, z_pxb, z_pxg, z_pxv)
    β_pxn = @. μ_pxn[u_n, :] + σ_pxn[u_n, :] * z_pxn[u_sn, :]

    # Likelihood
    β = sum(β_pxn[idx_sn, :] .* PC, dims=2)
    y ~ arraydist(LazyArray(@~ BernoulliLogit.(α_a[a] + α_s[s] + α_r[r] + α_sxr[idx_sr] + β)))

    # Generated quantities
    return (; y, α_a, α_s, α_r, α_sxr, β_pxn)
end;

# Same inputs as main model allows us to pass the same vector
function simulate(params, a, r, s, n, X, idx_sn, u_n, u_sn, Nv, Ng, Nb; idx_sr=idx(s, r))
    map(params) do p
        β = sum(p.β_pxn[idx_sn, :] .* X, dims=2)
        @. logistic(p.α_a[a] + p.α_s[s] + p.α_r[r] + p.α_sxr[idx_sr] + β)
    end
end

end

module CountModels

using Turing, LazyArrays, StatsFuns, LinearAlgebra
include("../src/utilities.jl")
import .CustomUtilityFuns: lu, idx

export mLc, simulate

@model function mLc(
    r, s, n, PC, y,
    idx_sn, u_n, u_sn, Nv, Ng, Nb;
    Nr=lu(r), Ns=lu(s), Nn=lu(n), NPC=size(PC, 2), idx_sr=idx(s, r)
)

    # Priors for species × region
    μ_sxr ~ filldist(Normal(0, σₚ), Ns)
    τ_sxr ~ filldist(InverseGamma(3, θₚ), Ns)
    z_sxr ~ filldist(Normal(), Ns, Nr)
    α_sxr = μ_sxr .+ τ_sxr .* z_sxr

    # Priors for nesting types × PCs
    μ_pxn ~ filldist(Normal(0, σₚ), Nn, NPC)
    τ_pxn ~ filldist(InverseGamma(3, θₚ/2), Nn, NPC)
    z_pxb ~ filldist(Normal(), Nb, NPC)
    z_pxg ~ filldist(Normal(), Ng, NPC)
    z_pxv ~ filldist(Normal(), Nv, NPC)
    z_pxn = ApplyArray(vcat, z_pxb, z_pxg, z_pxv)
    β_pxn = μ_pxn[u_n, :] .+ τ_pxn[u_n, :] .* z_pxn[u_sn, :]

    # Prior for random error
    σ2 ~ InverseGamma(3, θₚ)

    # Likelihood
    μ = vec(α_sxr[idx_sr] + sum(β_pxn[idx_sn, :] .* PC, dims=2))
    Σ = σ2 * I
    y ~ MvNormal(μ, Σ)

    # Generated quantities
    return (; y, α_sxr, β_pxn)
end;

function simulate(α, β, σ2, idx_sn, s, r, X; idx_sr=idx(s, r))
    out = Vector(undef, length(α))
    for i in eachindex(α)
        μ = α[i][idx_sr] .+ sum(β[i][idx_sn, :] .* X, dims=2)
        out[i] = rand.(Normal.(μ, σ2[i]))
    end
    return out
end
    
end