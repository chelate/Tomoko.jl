##

using Revise
##
using Tomoko

##

β1 = [ (mod(ii, 10)==0 ? 0.01 : 0) for ii in 1:256]
variable_sites = [ ii for ii in 1:256 if (mod(ii, 20)==0)]
flip_prob = .01;

##
par = PopRates(
    κ= 5000,
    χ=.2, 
    ρ=0.1, 
    μ=.1/5000, 
    loci = 256, 
    β1 = β1)
    # β1 is symmetrical 
pop = initialize_pop(par);
##
dg = run_sim(pop, par, 0:1:200 ; var_sites = variable_sites, flip_prob = 0)


##

using Plots
##
plot(dg.time, dg.mean_fit)