export PopRates, PopState, FitVector
export initialize_pop, random_wright_sample

## difinition of the fitness structure

struct FitVector # Custom sparse array type for fitness vector stuff
    β0::Float64
    sind::Array{Int64,1} # sites under selection
    nind::Array{Int64,1} # neutral sites
    svals::Array{Float64,1} # Δ fitness if mutant
    βwt::Float64 # β0 + sum(β1), the wildtype fitness (zero bitvector)
    βmax::Float64 # β0 + sum(abs.(β1)) the maximum fitness
end

function FitVector(β0::Float64, β1::Array)
    # vector of fitness differences between wildtype and mutant
    FitVector(
        β0,
        findall(!iszero, β1),
        findall(iszero, β1),
        -β1[findall(!iszero, β1)],
        β0 + sum(β1)/2,
        β0 + sum(abs.(β1))/2
    )
end


function fitness(f::FitVector, x::BitArray)
    s = 0
    for (ii,val) in zip(f.sind,f.svals)
        if x[ii]
            s += val
        end
    end
    return s + f.βwt
end


# the fitness vector is stored in cachey sparsy form.
# unfortunately we may want to access it in a more conventional way
function β1_vector(f::FitVector)
    # vector of fitness differences f(wildtype) - f(mutant)
    a = zeros(Float64, maximum([maximum(f.sind),maximum(f.nind)]))
    for (ii,val) in zip(f.sind,f.svals)
        a[ii] = -val
    end
    return(a)
end



##
# Paramters defining simulation state

mutable struct PopState{T}        # parameters defining a population state
    individuals::Vector{BitArray{1}}
    fitness::Vector{T}         # tendency toward excess birth
    birth::Vector{T}
    death::T                 # value shared among individuals
    size::Int
    birth_total::T            
    time::T
end

# Parameters defining the simulaiton dynamics
struct PopRates        # parameters defining a population
    loci::Int64        #    = 2^8        number of loci on the genome
    κ::Float64         #    = 1000.        avg pop size
    σ::Float64        #    = 0.5         excess Poisson noise
    μ::Float64         #    = 1.0e-3    mutation rate per site
    χ::Float64        #    = 0.0        outcrossing rate (between 0 and one), rate of recombination events
    ρ::Float64         #    = 1.0e-2    how tightly two crossed genomes get wound together (# of crosses per nucleotide)
    β0::Float64     #    = 1.0        base (competetive) fitness      
    βf::FitVector    #    = FitVector(β0,zeros(loci)) encodes fitness function efficiently
    f::Function     #    = (x->fitness(βf,x))    # mapping from genotype to fitness (increased birth + competion)
    # if you want epistasis, you can always code it in by hand in this function.
    function PopRates(loci,κ,σ,μ,χ,ρ,β0,β1)
        # we make βf and f only accessible through inner constructors
        # this enforces their correctness
        if length(β1) == loci
            βf = FitVector(β0,β1)
            f = (x->fitness(βf,x))
            return new(loci,κ,σ,μ,χ,ρ,β0,βf,f)
        else
            error("elSection coefficients wrong length")
        end
    end
end

# default constructor method
function PopRates(;
    loci = 2^8, κ = 1000, σ = 0.5, μ = 1.0e-3, χ = 0.0, ρ = 1.0e-2, β0 = 1.0,
    β1 = zeros(loci))
    PopRates(loci,κ,σ,μ,χ,ρ,β0,β1)
end

# updating method
function PopRates(par::PopRates; 
    loci = par.loci, κ = par.κ, σ = par.σ, μ = par.μ, χ = par.χ, ρ = par.ρ, β0 = par.β0,
    β1 = β1_vector(par.βf))
    PopRates(loci,κ,σ,μ,χ,ρ,β0,β1)
end


## Equilibrium Sampling
# We assume that the population is very near the maximum


function pop_ne(par::PopRates)
    par.κ/(2*(par.σ + par.βf.βmax))
end

function pop_D(par::PopRates)
    4*pop_ne(par)*par.μ*(par.σ + par.βf.βmax)
end

"""
Sample from Wright equilibrium
"""
function random_wright_sample(d, q, s)
    k = 0 
    while true
        k = rand(Poisson(s))
        p = (gamma(d)*gamma(k + d*q)) / (gamma(k + d)*gamma(d*q))
        rand(Bernoulli(p)) ? break : nothing
    end
    rand(Beta(q*d + k, (1 - q)*d)) # frequency
end

function wright_sample(par::PopRates)
    freq = zeros(par.loci)
    d = pop_D(par)
    ne = pop_ne(par)
    # set the selected frequencies
    for (ii, val) in zip(par.βf.sind, par.βf.svals)
        if val < 0 # if mutant has a disadvantage
            freq[ii] = 1 -random_wright_sample(d, 1//2, 2*ne*abs(val))
        else # if the mutant has an advantage
            freq[ii] = random_wright_sample(d, 1//2, 2*ne*abs(val))
        end
    end
    for ii in par.βf.nind
        freq[ii] = random_wright_sample(d, 1//2, 0)
    end
    return freq
end

"""
Construct a population individuals with site entries equal to 1 according to binomial
frequencies 
"""
function initialize_pop(frequencies::Vector, pop_size::Int, par::PopRates; time = 0.0)
    genomes = mapreduce(p->rand(Bernoulli(p), pop_size),hcat,frequencies)
    individuals = [BitArray(genomes[ii,:]) for ii in 1:pop_size]
    fitness = par.f.(individuals)
    birth = fitness .+ par.σ
    birth_total = sum(birth)
    death =  par.σ .+ sum(fitness)/par.κ
    pop = PopState(
        individuals, 
        fitness,
        birth,
        death,
        pop_size,
        birth_total,
        time)
end

# asumption of neutrality
function initialize_pop(par::PopRates; pop_size = ceil(Int, par.κ))
    initialize_pop(wright_sample(par), pop_size, par) # draw each site from Wright eq distribution
end

# simpler specification for constant frequencies
# useful for indicating pure wildtype
function initialize_pop(freq::Real, par::PopRates; pop_size = ceil(Int, par.κ))
initialize_pop([freq for ii in 1:par.loci], pop_size, par)
end

"""
Recalculate the caches
"""
function renew_fitness(pop::PopState, par::PopRates)
    pop.fitness = par.f.(pop.individuals)
    pop.birth = pop.fitness .+ par.σ
    pop.death = sum(pop.fitness)/par.κ .+ par.σ
    pop.size = length(pop.individuals)
    pop.birth_total = sum(pop.birth)
end

