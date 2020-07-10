##

using Tomoko
using Gadfly
##
par = PopRates(
    κ= 5000,
    χ=.2, 
    ρ=0.1, 
    μ=.1/5000, 
    loci = 256, 
    β1 = [ (mod(ii, 10)==0 ? 0.1 : 0) for ii in 1:256])
    # β1 is symmetrical 
pop = initialize_pop(0.0,par)
dg = run_sim(pop, par, 0:5:2000)
@profile df = run_sim(pop, par, 0:5:2000)