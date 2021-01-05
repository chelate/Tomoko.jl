using Distributed
addprocs(8)
@everywhere using Pkg
@everywhere Pkg.activate(".")
@everywhere using Tomoko
using CSV, DataFrames, Statistics


path = "/Users/lamont/Dropbox/Colin_ControlTheory/HIV trapping code/Julia results/"

## Check the sampler
vals = [Tomoko.sample!(Tomoko.WrightSampler(0.05,0.05,k*10)) for _ in 1:10^6, k in 0:6]
CSV.write(string(path,"sample",".csv"),DataFrame(vals),writeheader=false)
est = [log(1-k_antibodies_extinction(.005,1,100000,[[1,2,exp(ii)*1000]],samples=10^6)) for ii in -6:6, _ in 1:1]
CSV.write(string(path,"estprob",".csv"),DataFrame(est),writeheader=false)



## add an escape propbability to the statistics of the data frame
@everywhere Tomoko.eval(:(
function list_log_ext(pop::PopState, par::PopRates)
    map( s -> empirical_log_ext(pop, par; sites = s),
    hcat([ # each sublist is a list of possible escape site indices for an antibody (these are the simple kind)
        [[[ii+1]],
        [[ii+1],[ii+2]],
        [[ii+1],[ii+2],[ii+3]]]
        for ii in 10:20:500]...) )
end
))
@everywhere push!(Tomoko.pop_stats,:list_log_ext)
## Linkage simulations of escape probability

loci = 2^9
θ_ = .05
mu = θ_ /1000 # the first number sets the diversity parameter θ = 2N_e
λ = 4
β1 = zeros(loci)
selection_vals =  collect( 10:20:500)

for ii in selection_vals
    β1[ii+1] += (λ * mu / 2) * ii # Multiples of the direct mutation rate
    β1[ii+2] += (λ * mu / 2) * ii
    β1[ii+3] += (λ * mu / 2) * ii
end
par = PopRates(
    κ= 1000,
    χ= 0.0, 
    ρ= 0.1, 
    μ= mu , # first number determines D
    loci = loci, 
    β1 = β1)

dglist = pmap(x->run_sim(par, 1000:20:2000) , 1:2000)

## Extract dglist frequency distribution
function collate!(counts, observations)
    max = length(counts)-1
    for ii in observations
        if ii < max
            counts[round(Int,ii)+1] += 1
        else 
            counts[max+1] += 1
        end
    end
end


let counts = [zeros(100) for ii in 1:512]
    for dg in dglist
        obs = hcat((dg.freq .* dg.pop_size)...)
        for ii in eachindex(counts)
            collate!(counts[ii], obs[ii,:])
        end
    end
    global countmat = hcat(counts...)
end

exportmat=vcat(β1',countmat)
CSV.write(string(path,"histogram",".csv"),DataFrame(exportmat),writeheader=false)

##


ext_mat = mean(map(y->mean(map(x->exp.(x),y)), map(x->x.list_log_ext, dglist)))

using Plots

sim_mat = [pmap(ii->
    k_antibodies_extinction( par.κ * par.μ/2  # = θ_ts / 2
        , kk
        , 6 * par.κ * par.β0 / par.λ  # => 2 Ne f =  2 * κ * β_0/ λ 
        , [[1,1,ii]]; samples = 10^6)
        , selection_vals) for kk in 1:3]
p = plot(log10.(selection_vals),log10.(1 .- permutedims(ext_mat,[2,1])))
plot!(p, log10.(selection_vals), log10.(1 .- hcat(sim_mat...)))

##  Almost analytic results  #############
   


# at physical paramaters 
pars1074 = [
        [0.7639963491553017, 0.38199817457765084, 505.9257358559],
        [2.5279926983106034, 0.2974109056836004, 302.00565652011676],
        [1.3079957406811853, 0.34121628017770056, 168.57885310136976]]

pars3BNC = [
    [1.0659990872888254, 1.0659990872888254, 125.715013338639], 
    [1.1979972618664763, 0.47919890474659055,222.1651518113006], 
    [0.5, 0.5, 126.83837392646784]
]

nevec = 10 .^ collect( -5:.25:0)
escapevec = hcat([pmap(ii->k_antibodies_extinction(ii, kk, 3*10^5, pars1074,samples=10^7), nevec)
    for kk in 1:5]...)
CSV.write(string(path,"10-1074",".csv"),DataFrame(hcat(nevec,escapevec)),writeheader=false)
##

##
nevec = 10 .^ collect( -5:.25:0)
escapevec = hcat([pmap(ii->k_antibodies_extinction(ii, kk, 3*10^4, pars3BNC,samples=10^7), nevec)
    for kk in 1:5]...)
CSV.write(string(path,"3BNC",".csv"),DataFrame(hcat(nevec,escapevec)),writeheader=false)

##
nevec = 10 .^ collect( -5:.25:0)
escapevec = hcat([pmap(ii->k_antibodies_extinction(ii, kk, 3*10^5, pars3BNC,samples=10^7), nevec)
    for kk in 1:5]...)
CSV.write(string(path,"3BNC_lowgrowth",".csv"),DataFrame(hcat(nevec,escapevec)),writeheader=false)

##

using Plots
plot(log10.(nevec),log10.(1 .- escapevec))

