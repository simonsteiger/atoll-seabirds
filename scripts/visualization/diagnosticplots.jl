module DiagnosticPlots

export plot_diagnostic

using DataFrames, Turing, StatsPlots

function plot_rhat(c::Chains; alpha=1)
    df = DataFrame(describe(c)[1])
    rhat = df.rhat
    params = string.(df.parameters)
    scatter(rhat, params, alpha=alpha)
end

function plot_diagnostic(c::Chains, var::String; alpha=1)
    # If Rhat should be plotted, call separate function
    var == "rhat" && return plot_rhat(c; alpha)

    # Wrap in DataFrame
    df = DataFrame(c)

    # Check if var is any of the covered measures, otherwise error
    var ∈ names(df) || throw("$(var) is not a diagnostic measure.")

    # Get variable
    diagnostics = df[:, var]
    # Determine chain length
    nsamples = length(diagnostics)
    # Make sequence for plotting
    seq = Int64.(1:nsamples*0.005:nsamples)
    # Calculate rolling mean
    μ = [mean(diagnostics[seq[i]:seq[i+1]]) for i in eachindex(seq[begin:end-1])]
    # Calculate rolling variance
    σ = [std(diagnostics[seq[i]:seq[i+1]]) for i in eachindex(seq[begin:end-1])]
    # Assemble plot
    plot(seq[begin:end-1], μ, lw=1.5, label="Rolling μ")
    plot!(seq[begin:end-1], μ .+ σ, fillrange=μ .- σ, fillalpha=0.2, alpha=0, c=1, label="Rolling σ")
    # Add a horizontal line at 0.65 if acceptance_rate is plotted
    var == "acceptance_rate" && hline!([0.65], c=:red, lw=1.5, ls=:dash, label=:none)
    # Label axes
    xlabel!("Iteration")
    ylabel!(titlecase(replace(var, "_" => " ")))
    ylims!(extrema(diagnostics))
end

end