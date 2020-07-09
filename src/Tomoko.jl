module Tomoko
##
using StatsBase
using Parameters
using Random
using Distributions
using SpecialFunctions
using JuliennedArrays
using DataFrames

include("genetics.jl")

export run_sim


# struct encoding current state of population

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

To undertand the code, here we have randomly chosen cuts at 2/5.
From positions 3-5 the codes are swapped as in this diagram.
--v---v---
--^---^---

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
Add an individual to the population.
Adjust the birth and death rates accordingly
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

"""
Delete an individual to the population.
Adjust the birth and death rates accordingly
"""
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


"""
Choose an individual from the population (daddy) to recombine with the individual at
mommy_index.
"""

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

"""
Measure diversity of neutral sites
	data is a matrix of counts ( 2 x n)
"""



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


##
end # module
