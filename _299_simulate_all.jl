# ==============================================================================
# MASTER SIMULATOR
# Executes separate data simulation scripts sequentially in isolated processes.
# ==============================================================================

scripts = [
    "_200_simulate_samples.jl",
    "_210_simulate_low_pop.jl",
    "_220_simulate_low_sample.jl"
]

println("=== Starting Sequential Execution of $(length(scripts)) Scripts ===\n")

for (i, script_name) in enumerate(scripts)
    if !isfile(script_name)
        println("ERROR: Could not find file '$script_name'")
        continue
    end

    println("[$i/$(length(scripts))] Executing: $script_name")
    println("---------------------------------------------------")
    
    # Measure time for execution
    elapsed = @elapsed begin
        # $(Base.julia_cmd()) ensures the usage of the exact same Julia executable 
        # that is running this master script.
        try
            run(`$(Base.julia_cmd()) $script_name`)
        catch e
            println("!!! Error running $script_name: $e")
        end
    end
    
    println("---------------------------------------------------")
    println("Finished $script_name in $(round(elapsed, digits=2)) seconds.\n")
end

println("=== All tasks completed. ===")