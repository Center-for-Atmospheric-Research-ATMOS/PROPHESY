# counting noise
"""
    countElectrons(S::Cdouble,σ::Cdouble=1.0)

    models the counting noise when counting the number of electron hitting the CCD sensor
    The noise is Poisson distributed
    S: the deterministic average number of electrons
    σ: is the variance of the electronique noise (e.g. dark current) modeled by a Gaussian distribution in the Poisson parameter
"""
function countElectrons(S::Cdouble,σ::Cdouble=1.0)
    μ_count = S + σ*randn()
    μ_count = max(μ_count,0.0);
    rand(Poisson(μ_count))
end
function countElectrons(S::Array{Cdouble,1},σ::Cdouble=1.0)
    μ_count = S + σ*randn(length(S))
    μ_count[μ_count.<0.0] .= 0.0
    rand.(Poisson.(μ_count))
end