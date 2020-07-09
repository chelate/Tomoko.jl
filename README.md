# Tomoko

Pop-gen simulator for Julia.  This is a forward-time, individual-based, bi-allelic haploid simulator suitable for testing high-dimensional models of fitness, linkage, and recombination. It is named in honor of Tomoko Ohta, pioneer of the neutral theory in molecular genetics. Major inspiration in the API is from FFpopsim.

Included are also several tools associated with Wright equilibrium including sampling, and estimation of population and site parameters.

# Parameters and running simulations
The dynamics of a population are determined by a [Parameters](https://github.com/mauro3/Parameters.jl) type `PopRates`.

`PopRates` has the following default definitions

```julia
@with_kw struct PopRates{R}		# parameters defining a population
	κ::R 				= 1000.		# avg pop size
	β0::R 				= 1.0		# base (competetive) fitness  
	σ::R 				= 0.5 		# excess Poisson noise
	μ::R 				= 1.0e-3	# mutation rate per site
    χ::R				= 0.0		# outcrossing probability (between 0 and one), 
	ρ::R 				= 1.0e-2	# how tightly two crossed genomes get wound together (# of crosses per nucleotide)
	f::Function 		= (x->0.0)	# mapping from genotype to fitness (increased birth + competion)
end
```

To run a simulation, you 
<<<<<<< HEAD
1. define an instance of a `PopRates` object, 
2. construct an initial population using a number of different constructors and then 
3. `run_sim(pop,par,timesteps)` timesteps tells you both how long to run the simulation and at which time-steps to gather statistics.
=======
    1. define an instance of a `PopRates` object, 
    1. construct an initial population using a number of different constructors and then 
    1. `run_sim(pop,par,timesteps)` timesteps at which to gather statistics.
>>>>>>> master



```julia
par = PopRates(χ=.2, ρ=0.1, μ=10^-4)
pop = initialize_pop(par,gene_size = 2048)
@profile df = run_sim(pop::Population, par::PopRates, 1:5:10000)
```

