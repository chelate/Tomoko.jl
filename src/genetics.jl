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

struct FitVector # Custom sparse array type for fitness vector stuff
	sind::Array{Int64,1} # sites under selection
	nind::Array{Int64,1} # neutral sites
	svals::Array{Float64,1} # fitness addition, if mutant
	βwt::Float64 # β0 + sum(β1), the wildtype fitness (zero bitvector)
	βmax::Float64 # β0 + sum(abs.(a)) the maximum fitness
end

function FitVector(a::Array, β0::Float64)
	# vector of fitness differences between wildtype and mutant
	FitVector(
		findall(!iszero, a),
		findall(iszero, a),
		-a[findall(!iszero, a)],
		β0 + sum(a)/2,
		β0 + sum(abs.(a))/2
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


##
@with_kw struct PopRates		# parameters defining a population
	loci::Int64				= 2^8
	κ::Float64 				= 1000.		# avg pop size
	β0::Float64 			= 1.0		# base (competetive) fitness  
	β1::Array{Float64,1}	= zeros(loci)	# fitness difference between sites
	βf::FitVector			= FitVector(β1,β0)
	σ::Float64				= 0.5 		# excess Poisson noise
	μ::Float64 				= 1.0e-3	# mutation rate per site
	χ::Float64				= 0.0		# outcrossing rate (between 0 and one), rate of recombination events
	ρ::Float64 				= 1.0e-2	# how tightly two crossed genomes get wound together (# of crosses per nucleotide)
	f::Function 			= (x->fitness(βf,x))	# mapping from genotype to fitness (increased birth + competion)
	# if you want epistasis, you can always code it in by hand in this function.
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
		vec = sum(sample(pop.individuals, k, replace=false))[sites]
		vec = permutedims(hcat(vec , k .- vec))
		neutral_D(vec)
end

"""
If fitnesses are known, we calculate diversity only from the neutral sites.
"""
function stat_D(pop::Population, par::PopRates; 
	k = length(pop.individuals)) # size of population to take counts one, defaults to all
	stat_D(pop, sites = par.βf.nind, k = k)
end

function stat_freq(pop::Population)
	mean(pop.individuals)
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
    fitness = par.f.(individuals)
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
function initialize_pop(par::PopRates; pop_size = ceil(Int, par.κ))
initialize_pop([neutral_wright_sample(par::PopRates) for ii in 1:par.loci], pop_size, par)
end

# simpler specification for constant frequencies
# useful for indicating pure wildtype
function initialize_pop(freq::Real, par::PopRates; pop_size = ceil(Int, par.κ))
initialize_pop([freq for ii in 1:par.loci], pop_size, par)
end

"""
Recalculate the caches
"""
function renew_fitness(pop::Population, par::PopRates)
	pop.fitness = par.f.(pop.individual)
	pop.birth = pop.fitness .+ par.σ
	pop.death = sum(pop.fitness)/par.κ .+ par.σ
	pop.size = length(pop.individual)
	pop.birth_total = sum(pop.birth)
end

