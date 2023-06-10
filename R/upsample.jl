using CSV, DataFrames, Chain, StatsBase

df_anst = @chain CSV.read("data/envscores.csv", DataFrame) begin
    subset(_, :species => x -> x .== "Anous_stolidus")
    select(_, Not([:cond, :region, :Column1]))
end

function balance!(data, draws)
    for i in draws
        append!(data, DataFrame(data[Int64(i), :]))
    end
end

function upsample(data, N)
    p = sum(data.presence)
    a = nrow(data) - p
    a_N = Int64(N - a)
    p_N = Int64(N - p)

    # Throw error if p_N or a_N below 1
    a_N < 1 || p_N < 1 && throw("N too low to upsample.")

    p_range, a_range = 1:p, 1:a
    p_draws, a_draws = StatsBase.sample(p_range, p_N), StatsBase.sample(a_range, a_N)

    p_df = data[data.presence.==1.0, :]
    a_df = data[data.presence.==0.0, :]

    for c in [[p_df, p_draws], [a_df, a_draws]]
        balance!(c[1], c[2])
    end

    return append!(p_df, a_df)
end
