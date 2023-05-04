
using Test
using LinearAlgebra
using XPSsampling

#=Functions to be tested 

    # sampling
    export samplePosterior, acceptSample, samplePosteriorMeanAndCov
    export samplePosteriorModelMargin, acceptSampleModelMargin # marginalization over the measurement operator space (or some approximation of it)
    export samplePosteriorMargin, acceptSampleMargin

=#

function test_corrCovariance()
    w = ones(Cdouble,10);
    cor_len = 3.0;
    Γprior = corrCovariance(w;cor_len=cor_len);
    (typeof(Γprior)==Matrix{Float64}) & (all(Γprior.>=0.0)) & (!isinf(sum(Γprior)))
end

function test_smoothnessCovariance()
    w = ones(Cdouble,10);
    cor_len = 3.0;
    Dsqrt, Γprior, Dprior = smoothnessCovariance(w;cor_len=cor_len);

    # conditions
    cond1 = (typeof(Γprior)==Matrix{Float64}) & (all(Γprior.>=0.0)) & (!isinf(sum(Γprior)))
    cond2 = (typeof(Dsqrt)==Matrix{Float64}) & (!isinf(sum(Dsqrt))) & (!isnan(sum(Dsqrt)))
    cond3 = (typeof(Dprior)==Matrix{Float64}) & (!isinf(sum(Dprior))) & (!isnan(sum(Dprior))) 

    # results
    cond1 & cond2 & cond3
end


@testset "Sampling matrices" begin
  @test test_corrCovariance()
  @test test_smoothnessCovariance()
end


# communication mechanism 

function test_comMec()
    Nr = 10
    w = ones(Cdouble,10);
    cor_len = 5.0;
    Γsqrt = real(sqrt(corrCovariance(w;cor_len=cor_len)));
    σB = [1.0; 1.0];
    psmooth = 0.99;
    x_curr = 10ones(Cdouble,Nr);

    # simulate a sample (better test: simulate many samples and check their statistics)
    x_prop = transmissionMechanism(x_curr,Γsqrt,σB;psmooth=psmooth)

    # results 
    (typeof(x_prop)==Vector{Float64}) & (all(x_prop.>=0.0)) & (!isinf(sum(x_prop)))
end


@testset "Communication mechanism" begin
    @test test_comMec()
  end