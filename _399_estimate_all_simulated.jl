# ==============================================================================
# MASTER RUNNER
# Runs estimations on all simulated data.
# ==============================================================================

# Configuration: Set the Intermediate Count for each script here.
# A random subset of trials will be run if Intermediate Count > 0.
# Set to 0 to run the full simulation.
# WARNING: ALL OUTPUT FILES IN THE TARGET FOLDERS WILL BE OVERWRITTEN IF THEY EXIST.
# THIS INCLUDES INTERMEDIATE OR FULL SIMULATION FILES DEPENDING ON WHICH IS RUN HERE.
# ==============================================================================
CONFIG = [
    # (Script Path, Intermediate Count)
    ("./_300_estimate_simulated.jl", 1),   
    ("./_310_estimate_low_pop.jl", 1),
    ("./_320_estimate_low_sample.jl", 1),
    ("./_330_estimate_beta.jl", 1)
]

# --- Execution Loop ---
for (script_path, count) in CONFIG
    println("\n" * "="^60)
    println("LAUNCHING: $script_path")
    println("Intermediate Count: $count")
    println("="^60)

    # Set the environment variable just for this specific command.
    # 'withenv' temporarily sets the ENV var, runs the code, then resets it.
    withenv("JULIA_INTERMEDIATE_COUNT" => string(count)) do
        # This spawns a new Julia process.
        # It is exactly the same as typing `julia -t 10 script_name.jl` in your terminal.
        run(`julia -t 10 $script_path`)
        
    end
    
    println("COMPLETED: $script_path")
end

println("\n" * "="^60)
println("ALL SIMULATIONS FINISHED")
println("="^60)
