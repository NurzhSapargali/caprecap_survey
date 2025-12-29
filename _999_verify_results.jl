"""
Script to verify that intermediate results match full results within acceptable tolerance.
"""

using CSV
using DataFrames

OUTPUT_FOLDER::String = "./_900_output/data/"

"""
    run_comparison(inter_csv::String, full_csv::String)

Compares intermediate results CSV with full results CSV.
Reports any divergences where N_hat differs by more than ±2.
"""
function run_comparison(inter_csv::String, full_csv::String)
    if !isfile(inter_csv) || !isfile(full_csv)
        println("File error: Either '$inter_csv' or '$full_csv' does not exist.")
        return
    end

    df_i = CSV.read(inter_csv, DataFrame)
    df_f = CSV.read(full_csv, DataFrame)
    
    # Join on unique identifiers
    m = innerjoin(
        df_i,
        df_f,
        on = [:trial, :T, :alpha, :N, :type],
        renamecols = "_inter" => "_full"
    )
    
    println("\n" * "="^80)
    println("VERIFICATION REPORT: $(basename(inter_csv))")
    println("="^80)
    
    # Logic: N_hat within +- 2
    divergences = filter(r -> abs(r.N_hat_inter - r.N_hat_full) > 2.0, m)

    if nrow(divergences) == 0
        println("SUCCESS: All $(nrow(m)) trials matched within ± 2.")
    else
        println("MATCHES: $(nrow(m) - nrow(divergences))")
        println("DIVERGENCES: $(nrow(divergences))")
        println("\nDIAGNOSTIC DATA FOR UNSTABLE ESTIMATORS AND TRIALS:")
        
        # Select relevant columns and show without additional package formatting
        report = select(
            divergences, 
            :type,
            :trial,
            :T,
            :N, 
            :a_hat_inter,
            :a_hat_full, 
            :N_hat_inter,
            :N_hat_full
        )
        println(report)
        println("\n" * "-"^80)
    end
end

simulated_folder = "$(OUTPUT_FOLDER)simulated/"
appendix_folder = "$(OUTPUT_FOLDER)appendix/"
# Add main simulation result file pairs here
data_file_pairs = [
    ("$(simulated_folder)estimates_0.5_intermediate.csv", 
     "$(simulated_folder)estimates_0.5.csv"), 
    ("$(simulated_folder)estimates_2.0_intermediate.csv", 
     "$(simulated_folder)estimates_2.0.csv")    
]

low_pop_folder = appendix_folder * "low_pop/"
low_sample_folder = appendix_folder * "low_sample/"
beta_bin_folder = appendix_folder * "beta_bin/"

# Add low population result file pairs here
append!(data_file_pairs, [
    ("$(low_pop_folder)estimates_0.5_low_pop_intermediate.csv", 
     "$(low_pop_folder)estimates_0.5_low_pop.csv"), 
    ("$(low_pop_folder)estimates_2.0_low_pop_intermediate.csv", 
     "$(low_pop_folder)estimates_2.0_low_pop.csv")    
])

# Add low sample result file pairs here
append!(data_file_pairs, [
    ("$(low_sample_folder)estimates_0.5_low_sample_intermediate.csv", 
     "$(low_sample_folder)estimates_0.5_low_sample.csv"), 
    ("$(low_sample_folder)estimates_2.0_low_sample_intermediate.csv", 
     "$(low_sample_folder)estimates_2.0_low_sample.csv")    
])

# Add beta-bin result file pairs here
append!(data_file_pairs, [
    ("$(beta_bin_folder)estimates_0.5_betabin_intermediate.csv", 
     "$(beta_bin_folder)estimates_0.5_betabin.csv"), 
    ("$(beta_bin_folder)estimates_2.0_betabin_intermediate.csv", 
     "$(beta_bin_folder)estimates_2.0_betabin.csv")    
])

# Run comparisons
for (inter_file, full_file) in data_file_pairs
    run_comparison(inter_file, full_file)
end
