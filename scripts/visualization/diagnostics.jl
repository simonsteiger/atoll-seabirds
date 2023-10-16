module DiagnosticPlots

export plot_acceptance_rate

using DataFrames, Turing, StatsPlots

function plot_acceptance_rate(c::Chains)
    acc_rate = DataFrame(c).acceptance_rate
    nsamples = length(acc_rate)
    seq = Int64.(1:nsamples*0.005:nsamples)
    μ = [mean(acc_rate[seq[i]:seq[i+1]]) for i in eachindex(seq[begin:end-1])]
    σ = [std(acc_rate[seq[i]:seq[i+1]]) for i in eachindex(seq[begin:end-1])]
    plot(seq[begin:end-1], μ, lw=1.5, label="Rolling μ")
    plot!(seq[begin:end-1], μ .+ σ, fillrange=μ .- σ, fillalpha=0.2, alpha=0, c=1, label="Rolling σ")
    xlabel!("Iteration")
    ylabel!("Acceptance rate")
    hline!([0.65], c=:red, lw=1.5, ls=:dash, label=:none)
    ylims!(0, 1)
end

function plot_rhat(c::Chains; alpha=1)
    df = DataFrame(describe(c)[1])
    rhat = df.rhat
    params = string.(df.parameters)
    scatter(rhat, params, alpha=alpha)
end

end