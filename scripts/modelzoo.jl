module ModelZoo
    
    using Turing, LazyArrays
    include("../src/utilities.jl")
    using .CustomUtilityFuns

    export m1p, m1p_predict, m2p, mXp, mXp_predict

    # Helper for model declaration
    lu(x) = length(unique(x))

    @model function m1p(r, s, y; Nr=lu(r), Ns=lu(s), idx_sr=idx(s, r))
        # Priors for species × region
        α_sxr ~ filldist(Beta(1, 2), Ns * Nr)
    
        # Likelihood
        y ~ arraydist(LazyArray(@~ Bernoulli.(α_sxr[idx_sr])))

        return (; α_sxr)
    end;

    function m1p_predict(α, s, r; idx_sr=idx(s, r))
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

    @model function mXp(
        r, s, n, PC, # inputs
        idx_sn, u_n, u_sn, Nv, Ng, Nb, # indexes, lengths, etc
        y; # outcome
        Nr=lu(r), Ns=lu(s), Nn=lu(n), NPC=size(PC, 2), idx_sr=idx(s, r)
    )
    
        # Priors for species × region
        α_sxr ~ filldist(Normal(0, 1), Ns * Nr)
    
        # Priors for nesting types × PCs
        μ_pxn ~ filldist(Normal(0, 1), Nn, NPC)
        σ_pxn ~ filldist(InverseGamma(3, 2), Nn, NPC)
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
    function mXp_predict(α, β, r, s, n, X, idx_sn, u_n, u_sn, Nv, Ng, Nb; idx_sr=idx(s, r))
        [rand.(BernoulliLogit.(α[i][idx_sr] .+ sum(β[i][idx_sn, :] .* X, dims=2))) for i in eachindex(α)]
    end

end