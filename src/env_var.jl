export selection_flip, run_sim

function selection_flip(par::PopRates; var_sites = [], prob = 0.1)
    # Generates a new parameter struct with variable sites flipped at random with probability = prob
    n = length(var_sites) # length of vectors with selection
    k = rand(Binomial(n,prob))
    inds = var_sites[sample(1:length(var_sites),k, replace=false)]
    β1 = β1_vector(par.βf)
    β1[inds] .*= -1
    return PopRates(par, β1 = β1)
end

function run_sim(pop::PopState, par::PopRates, timepoints; var_sites = [], flip_prob = 0.1)
    df = initialize_stats()
    for t in timepoints
        run_until!(pop,par,t)
        record_stats!(df, pop, par)
        par = selection_flip(par; var_sites = var_sites, prob = flip_prob)
    end
    return df
end


function run_sim(par::PopRates, timepoints; var_sites = [], flip_prob = 0.1)
    pop = initialize_pop(par) # population loci set to wright equilibrium, pop.time set to zero
    df = initialize_stats() # empty data frame for holding the stats
    for t in timepoints
        run_until!(pop,par,t) # evolve the poulation state until pop.time > t.
        record_stats!(df, pop, par)  # record the statistics in the data_frame
        par = selection_flip(par; var_sites = var_sites, prob = flip_prob) # possibly stochastically change the environment.
    end
    return df
end


