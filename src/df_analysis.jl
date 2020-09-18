export estimate_Δ

"""
This is a MLE for fitting the selection parameters to a time trace, given the neutral parameters
"""

##
"""
The log-likelihood difference is computed using multiple importance sampling 
with a mixture of the neutral-beta and the binomial-reweighted neutral sample
1/2 * B(Dq, D(1-q)) + 1/2 * B(Dq+k, D(1-q) + N-k).
This give the benefit of using the same samples for both likelihoods, greatly reducing the noise.

"""

logbinomial(n, k) =  loggamma(1+n-k) - loggamma(1+n) - loggamma(1+k)
logbeta(α, β) = loggamma(α) + loggamma(β) - loggamma(α+β)

##
function randlogGamma(a,n;scale = 1)
    if a > .2
        return log.(rand(Gamma(a),n))
    else
        L = 1 / a - 1
        w = a / MathConstants.e / (1 - a)
        ww = 1 / (1 + w)
        η = z -> ( z >= 0 ? exp(-z) : w * L * exp(L * z))
        h = z -> exp(-z - exp(-z / a))
        function rh(a)
            while true
                U = rand()
                z = (U <= ww) ? -log(U / ww) : log(rand()) / L
                h(z) / η(z) > rand() && return(z)
            end
        end
        Z = Array{Float64}(undef,n)
        for ii in eachindex(Z)
            Z[ii] = log(scale) - rh(a) / a
        end 
        return Z
    end
end


function Δ_likelihood(D,q,n,k)
    n_samp = 1000
    loggam_a = vcat(
        randlogGamma(D*q, n_samp),
        randlogGamma(D*q+k, n_samp)
    )
    loggam_b = vcat(
        randlogGamma(D*(1-q), n_samp),
        randlogGamma(D*(1-q)+(n-k), n_samp)
    )
    xs = 1 ./(1 .+ exp.(loggam_b .- loggam_a))
    δ_z = logbeta(D*q, D*(1-q)) - logbeta(D*q+k, D*(1-q)+(n-k))
    δ_ls = map((a,b) -> k*a + (n-k)*b - n*(max(a,b) + log1p(exp(-abs(a-b)))) + δ_z, loggam_a, loggam_b) # log(p(d+k)/p(d)) 
    Δ -> log(sum( exp(x*Δ*D/2) / (1 + exp(δ_l)) for (x,δ_l) in zip(xs,δ_ls))) -
        log(sum( exp(x*Δ*D/2) / (1 + exp(-δ_l)) for (x,δ_l) in zip(xs,δ_ls))) - 
        logbinomial(n,k)
end


function make_likelihood(df; locus = 1, timepoints = 100:100:1000, sub_n = 100)
    ind = map( t->findfirst(x -> x>t, df.time),timepoints)
    ind = ind[ind.!=nothing]
    if length(ind) <= 1
        return (x->0)
    end
    Ds = df.D[ind] # population diversity
    # subsampling element (this may be neccesary to limit the large number of genomes)
    Ns = df.pop_size[ind] # total number of indivudals
    Ks = round.(Int,map(x -> x[locus], df.freq[ind]) .* Ns) #number of mutants.
    Js = Ns .- Ks
    ks = map((k,j) -> rand(Hypergeometric(k, j, min(sub_n,k+j))), Ks, Js)
    # if pop.size is smaller than the subsample size (sub_n), than default to pop.size
    fn_list = map((k,d,n)->Δ_likelihood(d,1//2,min(sub_n,n),k),ks,Ds,Ns)
    fn_normalize = Δ_likelihood(10^-5,1//2,2,1)
    # this function creates a minimum at finite Δ even if there is no mutants observed
    # D sets rough scale at which fitness effect will be invisible.
    Δ -> sum(f(Δ) for f in fn_list) + 0.01*fn_normalize(Δ)
end

function estimate_Δ(dfs::Vector; locus = 1, timepoints = 100:100:1000, sub_n = 100)
    fnlist = map(x-> make_likelihood(x; locus = locus, timepoints = timepoints, sub_n = sub_n),dfs)
    like(Δ) = mapreduce(f->f(Δ...), +, fnlist)
    s0 = [0.0]
    y = Optim.minimizer(optimize(like, s0, LBFGS(), autodiff = :forward))
    return -y[1]
end