# Tomoko

Pop-gen simulator for Julia.  This is a forward-time, individual-based, bi-allelic haploid simulator suitable for testing high-dimensional models of fitness, linkage, and recombination. It is named in honor of Tomoko Ohta, pioneer of the nearly-neutral theory in molecular genetics.  Major inspiration for the API is from the `high_d` simulator in [FFPopSim](http://webdav.tuebingen.mpg.de/ffpopsim/). 

Also included are several tools associated with Wright equilibrium including sampling and estimation of population and site parameters associated with the equilibrium distribution.

>Genetics is now at a very interesting stage. There are so many interesting questions unanswered and so many ways to test to find answers. Intuition is very important in addressing questions. Nurture your own sensibility and pursue your research and work with confidence. - Tomoko Ohta

The goal of this package is to help peope nuture their own sensibility about genetics and pursue their work with confidence.

# Design and limitations 

The code is very simple. Genotypes and individuals in a 1-1 mapping. We do not keep track of unique genotypes. Instead, we avoid O(`pop.size^2`) time scaling by using a rejection sampling scheme. The rejection sampling is most efficient when the nearly-neutral assumption holds, and the majority of the population has additive fitness close to the maximum.

# Parameters and running simulations
The dynamics of a population are determined by a [Parameters](https://github.com/mauro3/Parameters.jl) type `PopRates`.

`PopRates` has the following default definitions

```julia
@with_kw struct PopRates		# parameters defining a population
	loci::Int64				= 2^8
	κ::Float64 				= 1000.		# avg pop size
	β0::Float64 			= 1.0		# base (competetive) fitness  
	β1::Array{Float64,1}	= zeros(loci)	# fitness difference between sites.
	σ::Float64				= 0.5 		# excess Poisson noise
	μ::Float64 				= 1.0e-3	# mutation rate per site
	χ::Float64				= 0.0		# outcrossing rate (between 0 and one), rate of recombination events
	ρ::Float64 				= 1.0e-2	# how tightly two crossed genomes get wound together (# of crosses per nucleotide)
	f::Function				# = fitness function using β1 
end
```

`.β1` is keyword which is passed an array of fitness differences between wildtype (locus = 0) and the mutant (locus = 1).  These fitnesses are calculated by in a symmertric (-1,1) way, by default, with no epistasis.  Epistasis can be defined by defining an arbitrary fitness function 

To run a simulation, you 
1. define an instance of a `PopRates` object, 
2. construct an initial population using a number of different constructors and then 
3. `run_sim(pop,par,timesteps)` timesteps tells you both how long to run the simulation and at which time-steps to gather statistics.

`run_sim' propagates the population forward in time using an exact Gillespie simulatior and collects statistics at the specified timepoints. It returns these statistics and time points as a dataframe.

```julia
par = PopRates(χ=.2, ρ=0.1, μ=10^-4)
pop = initialize_pop(par)
df = run_sim(pop, par, 1:5:10000)
```
The statistics are scalar functions of the population. They are passed to `run_sim` as an optional vector of functions `stats = [f1, f2, ...]` and their corresponding symbols `names = [:f1,:f2,...]` then define the data-frame column names. 

The stats can be used for plotting purposes. For example,plotting the mean fitness over time from the previous example:
```julia
using Gadfly
plot(x = df.t, y = df.fit)
```
