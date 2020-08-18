


using Distributed
addprocs(7)
@everywhere using Pkg
@everywhere Pkg.activate(".")
@everywhere using Tomoko
using CSV
using Printf
using DataFrames
##
##

function save_pop_fit(path, result, par; name_fields = [:κ,:μ,:χ], p = 0)
    s = ["$(name)=$(@sprintf("%.0E", getfield(par,name)))_" for name in name_fields]
    CSV.write(string(path,s...,"pflip=$(@sprintf("%.0E", p))",".csv"),result)
end

path = "/Users/lamont/Dropbox/Colin_ControlTheory/HIV trapping code/Julia results/"


## One site under selection
loci = 2^10
# β1 = zeros(loci)
# fixed_sites = [ ii for ii in 1:loci if (mod(ii, 5)==0)]
# β1[fixed_sites] .= 0.00001 .* fixed_sites
# variable_sites = [ ii for ii in 1:loci if (mod(ii, 20)==3)]
# β1[variable_sites] .= 0.00001*5000
β1 = zeros(loci)
β1[50] += 5.0e-6 * 100
par = PopRates(
    κ= 1000,
    χ=0, 
    ρ=0.1, 
    μ= 0.01 /1000/2, # the first number sets the diversity 
    loci = loci, 
    β1 = β1)
dglist = pmap(x->run_sim(par, 0:10:5000 ; flip_prob =0),1:200)
##
fitness = pmap(ii->estimate_Δ(dglist, locus = ii, timepoints = 1000:100:5000), 40:60)
##



## Reference Simulation
global p = []
for D in [0.01,0.1]
par = PopRates(
    κ= 1000,
    χ=0, 
    ρ=0.1, 
    μ= D /1000/2, # the first number sets the diversity 
    loci = loci, 
    β1 = β1)
    # β1 is symmetrical 
true_fit =β1[fixed_sites]/par.μ
dglist = pmap(x->run_sim(par, 0:10:5000 ; flip_prob =0),1:ceil(Int64,(10/sqrt(D))))
fitness = pmap(ii->estimate_Δ(dglist, locus = ii, timepoints = 1000:100:5000), fixed_sites)
push!(p,[true_fit,fitness])
end

global q = []
for D in [0.01,0.1]
par = PopRates(
    κ= 1000,
    χ=1, 
    ρ=0.1, 
    μ= D /1000/2, # the first number sets the diversity 
    loci = loci, 
    β1 = β1)
    # β1 is symmetrical 
true_fit =β1[fixed_sites]/par.μ
dglist = pmap(x->run_sim_total_recombine(par, 0:10:5000 ; flip_prob =0),1:ceil(Int64,(10/sqrt(D))))
fitness = pmap(ii -> estimate_Δ(dglist, locus = ii, timepoints = 1000:100:5000),fixed_sites)
push!(q,[true_fit,fitness])
end


##
using Gadfly
ii = 1
plot(layer(x = p[ii][1], y = p[ii][2],Geom.line),
layer(x = q[ii][1], y = q[ii][2], Geom.line),
layer(x -> x, 0, maximum(p[ii][1])), Coord.cartesian(xmax =  maximum(p[ii][1]), ymax = 1.1*maximum(p[ii][2])))
##
par = PopRates(
    κ= 5000,
    χ=0.0, 
    ρ=0.1, 
    μ=.1 /5000/2, # the first number sets the diversity 
    loci = loci, 
    β1 = β1)
    # β1 is symmetrical 
true_fit =β1[fixed_sites]/par.μ
##

for x in [0,0.1], p in [0,0.01]
    result = DataFrame(f = true_fit)
    par = PopRates(
        κ= 5000,
        χ= x, 
        ρ= 0.1, 
        μ= .1 /5000/2 , # first number determines D
        loci = loci, 
        β1 = β1)
    
    for ii in 1:64
        dglist = [run_sim(par, 0:10:2000 ; var_sites = variable_sites, flip_prob = p/length(variable_sites)) for ii in 1:12]
        fitness = [estimate_Δ(dglist, locus = ii, timepoints = 1000:100:2000) for ii in fixed_sites]
        result[Symbol("run",ii)]=fitness
    end
    save_pop_fit(path,result,par, p = p)
end


##
df = run_sim(par, 0:20:13000 ; var_sites = variable_sites, flip_prob = .01/length(variable_sites))
##

l1 = layer(x =df.time ,y = (df.mean_fit .- mean(df.mean_fit) )./ 2 .+ 0.1 , Geom.line)
l2 = layer(x =df.time ,y = df.D , Geom.line, Theme(default_color=colorant"orange"))

plot(l1,l2, Coord.cartesian(xmax = 1e4))

##
dglist = [run_sim(par, 0:20:10000 ; var_sites = variable_sites, flip_prob = .01/length(variable_sites)) for ii in 1:4]
##
estimate_Δ(dglist, locus = 3, timepoints = 100:100:1000)

##

