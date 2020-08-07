
using Tomoko
##
##

β1 = [ (mod(ii, 10)==0 ? 0.0001*ii : 0) for ii in 1:256]
variable_sites = [ ii for ii in 1:256 if (mod(ii, 20)==0)]
flip_prob = .01;

##
par = PopRates(
    κ= 5000,
    χ=0.0, 
    ρ=0.1, 
    μ=0.05/5000, 
    loci = 256, 
    β1 = β1)
    # β1 is symmetrical 
pop = initialize_pop(par);
##
dglist = [run_sim(par, 0:20:1000 ; var_sites = variable_sites, flip_prob = 0) for ii in 1:4]
##
estimate_Δ(dglist, locus = 3, timepoints = 100:100:1000)

##

function estimate_Δ1(dfs; kwdargs...)
    fnlist = map(x-> Tomoko.make_likelihood(x,kwdargs...),dfs)
    like(Δ) = mapreduce(f->f(Δ...), +, fnlist)
    s0 = [0.0]
    y = Optim.minimizer(optimize(like, s0, LBFGS()))
    return y[1]
end
##
using Plots
##
pop = initialize_pop(par);
dg1 = run_sim(pop, par, 0:5:2000 ; var_sites = variable_sites, flip_prob = 0.01)
plot(dg1.time, dg1.mean_fit)
##
pop = initialize_pop(par);
dg2 = run_sim(pop, par, 0:5:2000 ; var_sites = variable_sites, flip_prob = 0.1)
plot(dg2.time, dg2.mean_fit)
pop = initialize_pop(par);
##
pop = initialize_pop(par);
dg3 = run_sim(pop, par, 0:5:2000 ; var_sites = variable_sites, flip_prob = 0.002)
plot(dg3.time, [0.2 .* dg3.mean_fit,dg3.D])


##
plot(dg1.time, dg1.D)
plot(dg2.time, dg2.D)
plot(dg3.time, dg3.D)
##

##
using HypergeometricFunctions

