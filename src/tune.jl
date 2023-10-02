using DataFrames

function tune(chain, test::Dict, test_label::Dict)
    df = DataFrame()

    for k in keys(test), i in 0:0.01:1
        predictions = prediction(test[k], chain, i)

        loss = sum((predictions - test_label[k]) .^ 2) / length(test_label[k])

        present = sum(test_label[k])
        absent = length(test_label[k]) - present

        predicted_present = sum(test_label[k] .== predictions .== 1)
        predicted_absent = sum(test_label[k] .== predictions .== 0)

        df = append!(
            df, DataFrame(
                "threshold" => i,
                "species" => k,
                "predicted_absent" => predicted_absent/absent,
                "predicted_present" => predicted_present/present,
                "criterion" => (predicted_absent/absent) * (predicted_present/present)
            ))
    end

    return df
end

function tune(chain, test, test_label)
    df = DataFrame()

    for i in 0:0.01:1
        predictions = prediction(test, chain, i)

        loss = sum((predictions - test_label) .^ 2) / length(test_label)

        present = sum(test_label)
        absent = length(test_label) - present

        predicted_present = sum(test_label .== predictions .== 1)
        predicted_absent = sum(test_label .== predictions .== 0)

        df = append!(
            df, DataFrame(
                "threshold" => i,
                "predicted_absent" => predicted_absent/absent,
                "predicted_present" => predicted_present/present,
                "criterion" => (predicted_absent/absent) * (predicted_present/present)
            ))
    end

    return df
end