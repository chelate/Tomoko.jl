

using Revise
using Distributed
addprocs(8)
@everywhere using Pkg
@everywhere Pkg.activate(".")
@everywhere using Tomoko
using Tomoko
using CSV
using Printf
using DataFrames
using Statistics


path = "/Users/lamont/Dropbox/Colin_ControlTheory/HIV trapping code/Julia results/"

function save_pop_fit(path, result, par; name_fields = [:κ,:μ,:χ], p = 0)
    s = ["$(name)=$(@sprintf("%.0E", getfield(par,name)))_" for name in name_fields]
    CSV.write(string(path,s...,"pflip=$(@sprintf("%.0E", p))",".csv"),result)
end

loci = 2^9
# β1 = zeros(loci)
# fixed_sites = [ ii for ii in 1:loci if (mod(ii, 5)==0)]
# β1[fixed_sites] .= 0.00001 .* fixed_sites
# variable_sites = [ ii for ii in 1:loci if (mod(ii, 20)==3)]
# β1[variable_sites] .= 0.00001*5000

##
for x in [0,0.1], p in [0,0.01], mu_ in [.1,.01]
    β1 = zeros(loci)
    mu = mu_ /1000/2 # the first number sets the diversity 
    fixed_sites = 2^3:2^3:2^9 # 64 fixed sites
    for ii in fixed_sites
        β1[ii] += mu * ii
    end
    variable_sites = 2^4-1:2^4:2^9 # 32 variable sites
    for ii in variable_sites
        β1[ii] += .02 # 2% fitness difference (max - mean < 4 to avoid blow up)
    end
    true_fit =β1[fixed_sites]/mu
    result = DataFrame(f = true_fit)
    par = PopRates(
        κ= 1000,
        χ= x, 
        ρ= 0.1, 
        μ= mu , # first number determines D
        loci = loci, 
        β1 = β1)
    
    for ii in 1:80
        dglist = pmap(x->run_sim(par, 0:10:2000 ; var_sites = variable_sites, flip_prob = p/length(variable_sites)), 1:12)
        fitness = pmap(ii->estimate_Δ(dglist, locus = ii, timepoints = 1000:100:2000) , fixed_sites)
        result[Symbol("run",ii)]=fitness
    end
    save_pop_fit(path,result,par, p = p)
end

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
using Gadfly
l1 = layer(x =df.time ,y = (df.mean_fit .- mean(df.mean_fit) )./ 2 .+ 0.1 , Geom.line)
l2 = layer(x =df.time ,y = df.D , Geom.line, Theme(default_color=colorant"orange"))

plot(l1,l2, Coord.cartesian(xmax = 1e4))

##
dglist = [run_sim(par, 0:20:10000 ; var_sites = variable_sites, flip_prob = .01/length(variable_sites)) for ii in 1:4]
##
estimate_Δ(dglist, locus = 3, timepoints = 100:100:1000)

##

