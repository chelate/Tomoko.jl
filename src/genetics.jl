export PopRates, Population
export initialize_pop, random_wright_sample


# Paramters defining simulation

mutable struct Population{T}		# parameters defining a population state
	individuals::Vector{BitArray{1}}
	fitness::Vector{T} 		# tendency toward excess birth
	birth::Vector{T}
	death::T 				# value shared among individuals
	size::Int
	birth_total::T			
end


@with_kw struct PopRates{R}		# parameters defining a population
	κ::R 				= 1000.		# avg pop size
	β0::R 				= 1.0		# base (competetive) fitness  
	σ::R 				= 0.5 		# excess Poisson noise
	μ::R 				= 1.0e-3	# mutation rate per site
	χ::R				= 0.0		# outcrossing rate (between 0 and one), rate of recombination events
	ρ::R 				= 1.0e-2	# how tightly two crossed genomes get wound together (# of crosses per nucleotide)
	f::Function 		= (x->0.0)	# mapping from genotype to fitness (increased birth + competion)
end


function PopRates(a::Array)
	PopRates(;f=(x->sum(a[ii]x[ii] for ii in eachindex(a))))
end

## Statistics on PopRates

function pop_ne(par::PopRates)
	par.κ/(2*(par.σ + par.β0))
end

function pop_D(par::PopRates)
	(4*par.κ*par.μ*(par.σ + par.β0))/(2*(par.σ + par.β0))
end

## Statistics on populations


function neutral_update_D(data::Array{T,2}, s::Real) where T
	# data has size = (number of categories, number of samples)
	mk = [0.5; 0.5] 
	ns = map(sum, Slices(data, 1))
	f1 = sum(digamma(s) - digamma(ns[ii] + s) +
			sum(mk[jj] * ( digamma(data[jj,ii] + s*mk[jj]) - digamma(s*mk[jj]))
			for jj in eachindex(mk))  
		for ii in eachindex(ns))
	f2 = sum(polygamma(1,s) - polygamma(1,ns[ii] + s) +
			sum(mk[jj]^2 * ( polygamma(1, data[jj,ii] + s*mk[jj]) - polygamma(1,s*mk[jj]))
				for jj in eachindex(mk)) for ii in eachindex(ns))
	a = -s^2*f2
	c = f1 - a/s
	if c >= 0
		a = s^3*(s*f2 + 2*f1)
		c = -(s^2 *f1 + a/s)
	end
	return -a/c
end

function neutral_D(data::Array{T,2}) where T
	s = 1
	while true
		s_new = neutral_update_D(data, s)
		abs(s - s_new) < 10^-9 ? break : s = s_new
	end
	return s
end

function stat_D(pop::Population; 
	sites = 1:length(pop.individuals[1]), # vector of sites to use in likelihood 
	k = length(pop.individuals)) # size of population to take counts one, defaults to all
		vec = sum(sample(pop.individuals, k, replace=false))
		vec = permutedims(hcat(vec , k .- vec))
		neutral_D(vec)
end

function stat_fit(pop::Population)
	mean(pop.fitness)
end

function stat_fitvar(pop::Population)
	var(pop.fitness)
end


##
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

function neutral_wright_sample(par::PopRates)
	random_wright_sample(pop_D(par), 1//2, 0)
end





"""
Construct a population individuals with site entries equal to 1 according to binomial
frequencies 
"""
function initialize_pop(frequencies::Vector, pop_size::Int, par::PopRates)
    genomes = mapreduce(p->rand(Bernoulli(p), pop_size),hcat,frequencies)
    individuals = [BitArray(genomes[ii,:]) for ii in 1:pop_size]
    fitness = par.β0 .+ par.f.(individuals)
    birth = fitness .+ par.σ
    birth_total = sum(birth)
    death =  par.σ .+ sum(fitness)/par.κ
    pop = Population(
        individuals, 
        fitness,
        birth,
        death,
        pop_size,
        birth_total)
end

# asumption of neutrality
function initialize_pop(par::PopRates; gene_size=128, pop_size = ceil(Int, par.κ))
initialize_pop([neutral_wright_sample(par::PopRates) for ii in 1:gene_size], pop_size, par)
end

# simpler specification for constant frequencies
function initialize_pop(freq::Real, par::PopRates; gene_size = 128, pop_size = ceil(Int, par.κ))
initialize_pop([freq for ii in 1:gene_size], pop_size, par)
end

"""
Recalculate the caches
"""
function renew_fitness(pop::Population, par::PopRates)
	pop.fitness = par.β0 .+ par.f.(pop.individual)
	pop.birth = pop.fitness .+ par.σ
	pop.death = sum(pop.fitness)/par.κ .+ par.σ
	pop.size = length(pop.individual)
	pop.birth_total = sum(pop.birth)
end

