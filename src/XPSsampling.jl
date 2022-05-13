# Sampling the a posteriori to estimate the uncertainty
#  - P(ρ|I,H): the probability density of the state ρ conditionally to the data I and the model H
#  - P(ρ|I):   the marginalization of P(ρ|I,H) over the measurement model space (using the small perturbation assumption)

# MAYBE: move to another package/module
# include("XPSsampling.jl")

module XPSsampling

    # covariance
    export smoothnessCovariance, corrCovariance

    # communication mechanism
    export transmissionMechanism

    # sampling
    export samplePosterior, acceptSample, samplePosteriorMeanAndCov
    export samplePosteriorModelMargin, acceptSampleModelMargin # marginalization over the measurement operator space (or some approximation of it)
    # export samplePosteriorEntropy, acceptSampleEntropy # does ot serve much purpose (except for showing that it's not gonna be used)
    export samplePosteriorMargin, acceptSampleMargin

    # implementation of the exported function/objects
    include("XPSsampling/sampleCov.jl")
    include("XPSsampling/communicationMechanism.jl")
    # include("XPSsampling/SAoptim.jl") # hide this for the time being
    # Metropolis Hastings samplers: fairly efficient since the sampled distributions are made up of the product of Gaussian distributions with a positivity constraint
    # [1] N. Metropolis et S. Ulam, The Monte Carlo method, 
    #     Journal of the American Statistical Association, vol. 44, no 247, 1949, p. 335–341 
    #     DOI: 10.2307/2280232
    # [2] W.K. Hastings, Monte Carlo Sampling Methods Using Markov Chains and Their Applications, 
    #     Biometrika, vol. 57, no 1, 1970, p. 97–109
    #     DOI: 10.2307/2334940
    # [3] wikipedia link https://en.wikipedia.org/wiki/Metropolis%E2%80%93Hastings_algorithm
    # include("XPSsampling/MHentropy.jl")      # hide this for the time being
    include("XPSsampling/MHquad.jl")         # good to use
    include("XPSsampling/MHquadBoundary.jl") # not really useful here

end


# """
#     Likelihood of the measurement knowing the state and the model

#     likelihood_H(I::Array{Cdouble,1},x::Array{Cdouble,1},H::Array{Cdouble,2},ΓIinv::Array{Cdouble,2},detΓI::Cdouble)

#     I: array of data
#     x: state of the system
#     H: measurement operator
#     ΓIinv: inverse of the measurement covariance matrix
#     detΓI: determinant of the measurement covariance matrix
# """
# function likelihood_H(I::Array{Cdouble,1},x::Array{Cdouble,1},H::Array{Cdouble,2},ΓIinv::Array{Cdouble,2},detΓI::Cdouble) # could become Poisson or something else
#     (1.0/sqrt(2π*detΓI))*exp(-0.5*(I-H*x)'*ΓIinv*(I-H*x))
# end

# """
#     accepting mechanism based on Metropolis-Hasting algorithm (for symmetric communication mechanisms)
# """

# """
#     Sampling the posterior distribution 
#     P(ρ|Y,H) ∝ P(Y|ρ,H)P(ρ)
#     retunr an array of samples ρ
# """