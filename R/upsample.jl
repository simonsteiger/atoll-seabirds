using CSV, DataFrames, Chain, DataFramesMeta, StatsBase

df_anst = @chain CSV.read("data/envscores.csv", DataFrame) begin
    subset(_, :species => x -> x .== "Anous_stolidus")
    select(_, Not([:cond, :region, :Column1]))
end

function balance!(data, draws)
    for i in draws
        append!(data, DataFrame(data[i, :]))
    end
end

function upsample!(data, N)
    p = sum(data.presence)
    a = nrow(data) - p
    N_a = Int64(N / 2 - a)
    N_p = Int64(N / 2 - p)

    # Throw error if N_p or N_a below 1
    N_a < 1 || N_p < 1 && throw("N too low to upsample.")

    p_range, a_range = 1:p, 1:a
    p_draws, a_draws = sample(p_range, N_p), sample(a_range, N_a)

    df_p = data[data.presence.==true, :]
    df_a = data[data.presence.==false, :]

    for c in [[df_p, p_draws], [df_a, a_draws]]
        balance!(c[1], c[2])
    end

    append!(df_p, df_a)
end

upsample(df_anst, 400)
