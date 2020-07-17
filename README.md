# Tomoko

Pop-gen simulator for Julia.  Tomoko is a forward-time, individual-based, bi-allelic haploid simulator suitable for testing high-dimensional models of fitness, linkage, and recombination. It is named in honor of Tomoko Ohta, pioneer of the nearly-neutral theory in molecular genetics.  Major inspiration for the API is from the `high_d` simulator in [FFPopSim](http://webdav.tuebingen.mpg.de/ffpopsim/). 

Also included are several tools associated with Wright equilibrium including sampling and estimation of population and site parameters associated with the equilibrium distribution.

Genetics is now almost 100 years old, but still very much active. Ohta herself has said,

>Genetics is now at a very interesting stage. There are so many interesting questions unanswered and so many ways to test to find answers. Intuition is very important in addressing questions. Nurture your own sensibility and pursue your research and work with confidence. 

This package was made to help nuture that sensibility.

# Parameters and running simulations
The scheme is simple. Genotypes (stored as `BitVector`'s) and individuals are in a 1-1 mapping. We do not keep track of the number of clones. Instead, we avoid O(`pop.size^2`) time scaling by using a rejection sampling scheme. The rejection sampling is most efficient when the nearly-neutral assumption holds, and the majority of the population has additive fitness close to the maximum.

The dynamics of a population are determined by a type `PopRates`.

`PopRates` has the following default definitions

```julia
struct PopRates		# parameters defining a population
	loci = 2^8      # Number of loci
	κ = 1000.       # avg pop size
	β0 = 1.0        # base (competetive) fitness  
	β1 = zeros(loci)# fitness difference between sites.
	σ = 0.5         # excess Poisson noise
	μ = 1.0e-3      # mutation rate per site
	χ = 0.0         # outcrossing rate (between 0 and one), rate of recombination events
	ρ  = 1.0e-2     # how tightly two crossed genomes get wound together (crosses per nucleotide)
	f::Function     # fitness function of genotype using sparse version of β0 + sum(β1 .*  x) 
end
```

`.β1` is keyword which is passed as an array of fitness differences between wildtype (locus = 0) and the mutant (locus = 1).  These fitnesses are implemented in a symmertric (-1,1) way by default, with no epistasis.  Epistasis can be defined by defining an arbitrary fitness function 

To run a simulation, you 
1. define an instance of a `PopRates` object, by specifying where the fields are different from the default values.
2. `run_sim(par,timesteps)` timesteps tells you both how long to run the simulation and at which time-steps to gather statistics.

```julia
par = PopRates(χ=.2, ρ=0.1, μ=10^-4, loci = 512)
df = run_sim(par, 1:5:10000)
```

`run_sim' initializes a population by drawing frequencies from the Wright equlibrium at each locus.  Then the population is propagated forward in time using an exact Gillespie simulatior while statistics at the specified timepoints. These statistics and time point are stored as a dataframe.

# The default statistics
We collect the following default statistics as columns in a DataFrame
`[:time, :pop_size, :freq, :D, :mean_fit, :var_fit]'

* `time` The time at which the statistics were collected. Interconvertable with generations in the Wright-Fisher sense
* `pop_size` Number of individuals at a given time
* `freq` vector of frequency of mutants (i.e. ones) at a given loci
* `D` The genotypic diversity of the population. (technically defined as the maximum likelihood estimator for the precision of the beta-binomial distribution of neutral site allellic counts.)
* `mean_fit` The mean fitness of the population
* `var_fit` The variance of the fitness of the population

Plotting, especially with `Gadfly.jl` and `StatPlots.jl` is well integrated with data frames.  For example, using the output of the simulation above, you can plot the `pop_size' over time

```julia
using Gadfly
plot(df, x = :time, y = :pop_size, Geom.line)
```

# Defining more complicated simulations.

In the simulation defined, we start with a populaiton in pseudo-equilibrium and observe the (mild) effects of linkage and drift when the simulation runs. A pseudo-equilibrium starting point is useful if we are interested in the population stationary state because we have to wait less time for things to settle down. However, you are probably running exact simulations because you are interested in what happens *far* from stationary state. This section is about how to spice things up a bit.

## What is inside `run_sim`

Let's just look at what `run_sim` is made of. It's nothing much more than a for loop. Here's the main definition

```julia
function run_sim(pop::Popstate, par::PopRates, timepoints;
	var_sites = [], # vector of indices that are flipping
	flip_prob = 0.1) # probability that they flip at every time point.
    df = initialize_stats() # empty data frame for holding the stats
    for t in timepoints
        run_until!(pop,par,t) # evolve the poulation state until pop.time > t.
        record_stats!(df, pop, par)  # record the statistics in the data_frame
        par = selection_flip(par; var_sites = var_sites, prob = flip_prob) # possibly stochastically change the environment.
    end
    return df # return the df, now heavy with juicy statistics for you to plot or analyze
end

# definition when pop is not provided: linkage-free equilibrium.
function run_sim(par::PopRates, timepoints; kwdargs...) # population loci set to wright equilibrium, 
	pop = initialize_pop(par)
	run_sim(pop::Popstate, par::PopRates, timepoints; kwdargs...)
end

```
Complications can fall into three categories
1. Change start population
1. Vary environment or add external events
1. Calculate more statistics.


## Shanging the start population
To make things more intersting, we can change the initial population to be more out-of equilibrium. For instance, we might define the population to be entirely made up of wildtype clones

```julia
pop = initialize_pop(0.0,par) # the first argument specifies the mutant frequency
```

Or we can grow out the population out from just a few individuals

```julia
pop = initialize_pop(0.0,par; pop_size = 10) 
par = PopRates(χ=.2, ρ=0.1, μ=10^-4)
df = run_sim(pop, par, 1:5:10000)
```
We can also define the population with a full vector of frequencies of length equal to the number of loci.

## Varying the evironment, bottleneck events.
As seen above, `run_sim` has the ability to include environmental variation, consisting of selection sign flips on some number of active sites through the defined `active_sites` vector.  This gets used by `selection_flip` which returns a new parameter value.

One can replace where `selection_flip` occurs in the `for` loop with more complicated changes to the evolutionary parameters or with functions on the population, like bottleneck events that remove individuals based on a particular sequence of loci.

Note: `par::PopRates` is an immutable with internal constructors and isn't be mutated on the fly. Instead a new instance is defined and overwrites the local variable `par`. On the other hand,`pop::PopState` is a mutable and can be mutated at will.

## Adding more statistics
The statistics are functions of the population and parameters that live inside the Tomoko module.  This is how `record_stats` identifies the statistic name with the statistic function.  This means that to define more functions from the REPL you have to eval them into the Tomoko context and update the global variable `pop_stats` so that `record_stats` knows you want to keep track of a new variable.

```
julia> Tomoko.eval(:(
function new_stat(pop::PopState, par::PopRates)
    do_stuff(pop::PopState, par::PopRates, sites = sites_of_interest)
end
))

julia> push!(Tomoko.pop_stats,:new_stat)
```

This design was chosen to keep the statistic-gathering machinery as global variables to avoid having to pass yet another struct to the simulation functions and to make sure that the name and function are intrinsincally linked.

# Extending and contributing
In the end it's impossible to design a user interface that can do everything from the REPl.   To run the experiment you need to run to answer your scientific questions, you will probably have to look at the source and see what's there, dev and modify the package to suit your needs.

By understanding what methods are available out of the box, you can get a feel for how the machinery works and how to extend it. Help Tomoko.jl evolve with PR's and feature requests!