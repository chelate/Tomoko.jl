using Tomoko
using Random
using Test

p = PopRates()
x = bitrand(p.loci)
df = runsim(p,1:1:1000)

@testset "Tomoko.jl" begin
    0.0 < p.f(x)
    (mean(df.D) - Tomoko.pop_D(p))/Tomoko.pop_D(p) < 0.05 # check that D is near the equilibrium value
    (mean(df.pop_size) - p.κ)/p.κ < 0.05 # check that the pop_size approaches the eq value
end