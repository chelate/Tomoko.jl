module Tomako
using StatsBase
using Parameters
using Random
using Distributions
using SpecialFunctions
using JuliennedArrays
using DataFrames

# Paramters defining simulation

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

# struct encoding current state of population

mutable struct Population{T}		# parameters defining a population state
	individuals::Vector{BitArray{1}}
	fitness::Vector{T} 		# tendency toward excess birth
	birth::Vector{T}
	death::T 				# value shared among individuals
	size::Int
	birth_total::T			
end

function  gillespie_step(pop::Population)
	total = (pop.birth_total + pop.death*pop.size)
	Δt = randexp() / total
	if rand() < (pop.birth_total / total)
		process = sample(1:pop.size, Weights(pop.birth, pop.birth_total))
	else
		process = - sample(1:pop.size)
	end
	return (process, Δt)
end

# Mutation

function flipat!(bits::BitArray, inds)
	for ii in inds
		bits[ii] = !bits[ii]
	end
end

"""
Mutates the genome at a rate given by μ
Assumes the action of a markov chain so that 
mutate!(X, 2μ) = mutate!(mutate!(X, μ), μ).
mutate!(Inf) should randomize the genome
"""

function mutate!(genome::BitArray, μ::Real)
	N = length(genome)
	p = -expm1(-2μ)/2
	flips = randsubseq(1:N, p)
	flipat!(genome, flips)
end

"""
Recombines two genomes at a rate given by μ
Assumes the action of a markov chain so that 
mutate!(X, 2μ) = mutate!(mutate!(X, μ), μ).
mutate!(Inf) should randomize the genome
"""


function recombine!(v1::T, v2::T, ρ::Real) where T <: AbstractArray
	N = length(v1)
	cut = sort(randsubseq(0:(N-1),ρ))
	if isodd(length(cut))
		push!(cut,N)
	end
	for ii in 1:2:length(cut)
	 	for jj in cut[ii]+1:cut[ii+1]
	 		(v2[jj], v1[jj]) = (v1[jj], v2[jj])  
	 	end
	end
end

"""
add an individual to the population, adjusting the birth and death rates accordingly
"""
function add_individual!(pop::Population, child::BitArray{1}, par::PopRates)
	# add to the fitness
	child_fitness = par.β0 + par.f(child)
	pop.birth_total += child_fitness + par.σ
	pop.death +=  child_fitness/par.κ			# add to the death rate
	push!(pop.birth, child_fitness + par.σ) 	# add to the birthrates
	push!(pop.fitness, child_fitness)
	push!(pop.individuals, child)
	pop.size += 1 
end

function delete_individual!(pop::Population, index::Int, par::PopRates)
	pop.birth_total += - pop.birth[index]
	# add to the death rate
	pop.death +=  -pop.fitness[index]/par.κ
	deleteat!(pop.birth, index)
	deleteat!(pop.fitness, index)
	deleteat!(pop.individuals, index)
	pop.size += -1 
end

function birth_event!(pop::Population, parent_index::Int, par::PopRates)
	child = copy(pop.individuals[parent_index])
	mutate!(child, par.μ)
	# add to the population
	add_individual!(pop, child, par)
end

death_event!(pop::Population, index::Int, par::PopRates) = delete_individual!(pop, index, par)

function recombination_event!(pop::Population, mommy_index::Int, par::PopRates)
	daddy_index = rand(1:pop.size)
	# update their genomes
	@views recombine!(pop.individuals[mommy_index], pop.individuals[daddy_index], par.ρ)
	# the fitness differences
	Δmfit = par.β0 + par.f(pop.individuals[mommy_index]) - pop.fitness[mommy_index] 
	Δdfit = par.β0 + par.f(pop.individuals[daddy_index]) - pop.fitness[daddy_index]
	pop.fitness[mommy_index] += Δmfit
	pop.fitness[daddy_index] += Δdfit
	# update birthrates
	pop.birth[mommy_index] += Δmfit
	pop.birth[daddy_index] += Δdfit
	pop.birth_total += Δmfit + Δdfit
	# update deathrate
	pop.death += (Δmfit+Δdfit)/par.κ
end

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

function next_event!(pop::Population, par::PopRates)
	# step forward in time
	(index, Δt) = gillespie_step(pop::Population)
	if index > 0
		birth_event!(pop, index, par)
		if rand(Bernoulli(par.ρ))
			recombination_event!(pop, index, par)
		end
	else
		death_event!(pop, abs(index), par)
	end
	Δt
end


## Genetic Functions

function pop_ne(par::PopRates)
	par.κ/(2*(par.σ + par.β0))
end

function pop_D(par::PopRates)
	(4*par.κ*par.μ*(par.σ + par.β0))/(2*(par.σ + par.β0))
end


## Function Goals
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
Measure diversity of neutral sites
	data is a matrix of counts ( 2 x n)
"""

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

function run_sim(pop::Population, par::PopRates, timepoints; 
	stats = [stat_D, stat_fit, stat_fitvar], 
	names = [:D, :fit, :fitvar] )
	rows = []
	time = 0 
	for t in timepoints
		while time < t
			time += next_event!(pop, par) 
		end
		push!(rows,[t;[s(pop) for s in stats]])
	end
	columns = Vector{Vector}(undef, 0)
	for (i,r) in enumerate(first(rows))
           column = Vector{typeof(r)}(undef, length(rows))
           column .= getindex.(rows, i)
           push!(columns, column)
       end
	df = DataFrame(columns, [:t;names], copycols=false)
end


end # module
