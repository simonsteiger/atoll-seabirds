using CSV, DataFrames, Chain, StatsBase

df_anst = @chain CSV.read("data/envscores.csv", DataFrame) begin
    subset(_, :species => x -> x .== "Anous_stolidus")
    select(_, Not([:cond, :region, :Column1]))
end

function balance!(data, draws)
    for i in draws
        append!(data, DataFrame(data[i, :]))
    end
end

function upsample(data, N)
    p = sum(data.presence)
    a = nrow(data) - p
    a_N = Int64(N / 2 - a)
    p_N = Int64(N / 2 - p)

    # Throw error if p_N or a_N below 1
    a_N < 1 || p_N < 1 && throw("N too low to upsample.")

    p_range, a_range = 1:p, 1:a
    p_draws, a_draws = sample(p_range, p_N), sample(a_range, a_N)

    p_df = data[data.presence.==true, :]
    a_df = data[data.presence.==false, :]

    for c in [[p_df, p_draws], [a_df, a_draws]]
        balance!(c[1], c[2])
    end

    return append!(p_df, a_df)
end

res = upsample(df_anst, 400)
