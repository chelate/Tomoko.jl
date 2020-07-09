# Tomoko

Pop-gen simulator for Julia.  This is a forward-time, individual-based, bi-allelic haploid simulator suitable for testing high-dimensional models of fitness, linkage, and recombination. It is named in honor of Tomoko Ohta, pioneer of the nearly-neutral theory in molecular genetics.  Major inspiration for the API is from the `highd` simulator in [FFPopSim](http://webdav.tuebingen.mpg.de/ffpopsim/). 

Included are also several tools associated with Wright equilibrium including sampling and estimation of population and site parameters associated with the equilibrium distribution.

# Parameters and running simulations
The dynamics of a population are determined by a [Parameters](https://github.com/mauro3/Parameters.jl) type `PopRates`.

`PopRates` has the following default definitions

```julia
@with_kw struct PopRates{R}	# parameters defining population dynamics
	κ::R 		= 1000.		# avg pop size
	β0::R 		= 1.0		# base (competetive) fitness 
	σ::R 		= 0.5 		# excess Poisson noise
	μ::R 		= 1.0e-3	# mutation rate per site
    χ::R		= 0.0		# outcrossing probability (0:1)
	ρ::R 		= 1.0e-2	# crossing per nucleotide
	f::Function = (x->0.0)	# mapping from genotype to fitness
end
```

To run a simulation, you 
1. define an instance of a `PopRates` object, 
2. construct an initial population using a number of different constructors and then 
3. `run_sim(pop,par,timesteps)` timesteps tells you both how long to run the simulation and at which time-steps to gather statistics.

`run_sim' propagates the population forward in time using an exact Gillespie simulatior and collects statistics at the specified timepoints.

```julia
par = PopRates(χ=.2, ρ=0.1, μ=10^-4)
pop = initialize_pop(par,gene_size = 2048)
df = run_sim(pop::Population, par::PopRates, 1:5:10000)
```
The statistics, passed to run_sim as an optional vector of functions, include the timestamp are stored in a data frame. Statistics are generally scalar functions of the population.

# Epilogue

Ohta tells us:

>Genetics is now at a very interesting stage. There are so many interesting questions unanswered and so many ways to test to find answers. Intuition is very important in addressing questions. Nurture your own sensibility and pursue your research and work with confidence.

Hopefully this package will help people test to find answers in Genetics.