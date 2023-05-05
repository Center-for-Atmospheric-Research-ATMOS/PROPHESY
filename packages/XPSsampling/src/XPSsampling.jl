# Sampling the a posteriori to estimate the uncertainty
#  - P(ρ|I,H): the probability density of the state ρ conditionally to the data I and the model H
#  - P(ρ|I):   the marginalization of P(ρ|I,H) over the measurement model space


"""
This is the [`XPSsampling`](@ref), which contains
* [`XPSsampling.smoothnessCovariance`](@ref)
* [`XPSsampling.corrCovariance`](@ref)
* [`XPSsampling.transmissionMechanism`](@ref)
* [`XPSsampling.samplePosterior`](@ref)
* [`XPSsampling.acceptSample`](@ref)
* [`XPSsampling.samplePosteriorMeanAndCov`](@ref)
* [`XPSsampling.samplePosteriorMargin`](@ref)
* [`XPSsampling.acceptSampleMargin`](@ref)
* [`XPSsampling.samplePosteriorBoundary`](@ref)
* [`XPSsampling.acceptSampleBoundary`](@ref)
* [`XPSsampling.samplePosteriorModelMargin`](@ref)
* [`XPSsampling.acceptSampleModelMargin`](@ref)
"""
module XPSsampling

    using LinearAlgebra

    # covariance
    export smoothnessCovariance, corrCovariance

    # communication mechanism
    export transmissionMechanism

    # sampling
    export samplePosterior, acceptSample, samplePosteriorMeanAndCov # MHquad.jl
    export samplePosteriorMargin, acceptSampleMargin
    export samplePosteriorBoundary, acceptSampleBoundary
    export samplePosteriorModelMargin, acceptSampleModelMargin      # MHquadBoundary.jl # approximation of the marginalization over the measurement operator space (low dimension)
    

    # low level functions: communication mechanisms and matrices
    include("sampleCov.jl")
    include("communicationMechanism.jl")
    # Metropolis Hastings samplers: fairly efficient since the sampled distributions are made up of the product of Gaussian distributions with a positivity constraint
    # [1] N. Metropolis et S. Ulam, The Monte Carlo method, 
    #     Journal of the American Statistical Association, vol. 44, no 247, 1949, p. 335–341 
    #     DOI: 10.2307/2280232
    # [2] W.K. Hastings, Monte Carlo Sampling Methods Using Markov Chains and Their Applications, 
    #     Biometrika, vol. 57, no 1, 1970, p. 97–109
    #     DOI: 10.2307/2334940
    # [3] wikipedia link https://en.wikipedia.org/wiki/Metropolis%E2%80%93Hastings_algorithm

    # sampling and estimations (covariance matrix)
    include("MHquad.jl")         
    include("MHquadBoundary.jl") # not really useful here

end
