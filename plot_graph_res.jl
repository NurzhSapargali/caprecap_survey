using DelimitedFiles
using DataFrames
using Plots

OUTPUT_FOLDER::String = "./_900_output/data/graphs/"
INPUT_FOLDER::String = "./_200_input/graphs/"
STRUCTURES::Vector{String} = ["nodes", "edges", "tris", "4cliques"]
EDGES_PER_NODE::Vector{Int64} = [1, 2, 3]

for en in EDGES_PER_NODE
    metafile = INPUT_FOLDER * "ba_$(en)/metadata.csv"
    for s in STRUCTURES
        file = OUTPUT_FOLDER * "ba_$(en)/$(s)_estimates.csv"
        data = try readdlm(file, ',', header=false) 
        catch error
            if isa(error, ArgumentError)
                continue
            end
        end

        df = DataFrame(data, ["alpha_hat", "Nu_hat", "No", "T", "avg_n"])
        df[!, "N_hat"] = df[!,"Nu_hat"] + df[!, "No"]
        units = "Nodes"
        if s == "edges"
            units = "edges"
        elseif s == "tris"
            units = "triangles"
        elseif s == "4cliques"
            units = "4-cliques"
        end
        p1 = histogram(df[!, "No"], title="$(titlecase(units)) discovered", label="")
        p2 = histogram(df[!, "Nu_hat"], title="Unobserved $(units) estimated", label="")
        p3 = histogram(df[!, "N_hat"], title="Total $(units) estimated", label="")
        vline!([sum(df[!, "N_hat"]) / length(df[!, "N_hat"])], label="Mean estimate", color="red", lw=5)
        meta = readdlm(metafile, ' ', header=false)
        truth = meta[1, findfirst(x -> x == s, STRUCTURES)]
        vline!([truth], label="True value", color="purple", lw=5)
        plot(p1, p2, p3, layout=(3, 1))
        estis = df[!, ["Nu_hat", "No", "N_hat"]]
        plot!(xlim=[minimum(minimum.(eachcol(estis))), 
                    maximum(maximum.(eachcol(estis)))])
        plot!(size=(1280, 720))
        savefig("$(s)_ba_$(en).pdf")
    end
end
