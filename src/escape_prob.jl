export k_antibodies_extinction, empirical_log_ext

# This sampler is robuset for small-exponent gamma distributions
# Code was copied from somewhere unknown (citation needed)
function randlogGamma(a;scale = 1)
    if a > .2
        return log.(rand(Gamma(a)))
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
        return log(scale) - rh(a) / a
    end
end

# This lists the neccesary parameters 
# for mcmc sampling from the Wright-equilibrium distribution
# αf, αb, σ are forward and backward diversity parameters respectively
mutable struct WrightSampler
    αf::Float64
    αb::Float64
    σ::Float64
    k::Int64
end

function WrightSampler(αf, αβ, σ)
    k = rand(Poisson(σ))
    ws = WrightSampler(αf, αβ, σ,k)
    for _ in 1:k 
    incrementk!(ws)
    end
    return ws
end

function incrementk!(ws::WrightSampler)
    state = ws.k
    step = rand([-1,1])
    next = state + step
    if next >= 0 #only make the change if in the positive domain#
        if step == 1
            p = ws.σ * (state + ws.αb) / (next * (state + ws.αf + ws.αb))
        elseif step == -1
            p = (state * (next + ws.αf + ws.αb)) / ws.σ * (next + ws.αb)
        end
        if rand(Bernoulli(min(p,1)))
            ws.k = next
        end
    end
end

function beta_draw(ws::WrightSampler)
    a = randlogGamma(ws.αb + ws.k;scale = 1)
    b = randlogGamma(ws.αf;scale = 1)
    1/(1+exp(b-a))
end

function sample!(ws::WrightSampler)
    incrementk!(ws)
    beta_draw(ws)
end

function escape_freq(freq_array)
    prod(1 .- prod(freq_array, dims = 1))
end

function prob_of_extinction(σf, sampler_array, samples)
    mean( exp(- σf * escape_freq(sample!.(sampler_array))) for ii in 1:samples)
end


# since ab_profile is measured in transitions, actually need ne to be 1/2 θ_ts
# θ_ts = 2 n/λ * (β *mu) = 2 n mu /2 = n mu
# fitness: σf = (2 * ne) .* f = θ_ts f  => f = β0/(β mu) =>
# This implies f = should be 2 β_0 /λ mu *
function k_antibodies_extinction(par::PopRates, k, ab_profile; samples = 10^5)
    k_antibodies_extinction( par.κ * par.μ / 2 # = θ_ts / 2
    , k
    , 2*par.β0 / (par.λ * par.μ) # => 2 Ne f =  2 * κ * β_0/ λ 
    , ab_profile; samples = samples)
end


# k is the number of antibodies
# f is aboslute growth rate (free fitness)
# ab_profile is fitness cost and mutational bias measured in transitions.
function k_antibodies_extinction(ne,k,f,ab_profile; samples = 10^5)
    sampler_array = map(x->WrightSampler(x...),
        hcat([(2 * ne) .* ab_profile for  j in 1:k]...))
    prob_of_extinction((2 * ne) .* f, sampler_array, samples)
end


function empirical_log_ext(pop::PopState, par::PopRates;
    sites = [[1,2],[3]] # must escape from and(or(1,2),3)
    )
    escape_fn(genome::BitArray) =  Bool(prod( 
        1 - prod(1-genome[ii] for ii in jj) # 1 if one of the sites is escaped
        for jj in sites)) # 1 if all the sites have atleast 1 escape
    survivors = filter(escape_fn, pop.individuals)
    birth_rates = (par.λ .+ par.f.(survivors)) ./ 2
    death_rates = (par.λ .- par.f.(survivors)) ./ 2
    log_ext = sum(log.(death_rates) .- log.(birth_rates))
end


# Only applicable to a single site, after simulation
# function df_log_ext(df, par::PopRates, fit;
#     sites = 1 # must escape from and(or(1,2),3)
#     )
#     out = log(1 - 2 fit / par.λ) * df.pop_size[] .* df.freq[ind]
# end