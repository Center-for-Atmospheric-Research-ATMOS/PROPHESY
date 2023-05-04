
using Test
using LinearAlgebra
using XPSsampling

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


function test_acceptSample()
    Nr = 10;
    ρ_cur  = 10ones(Cdouble,Nr);
    ρ_prop = 9ones(Cdouble,Nr);
    H = diagm(0 => ones(Cdouble,Nr) ,1 => 0.5ones(Cdouble,Nr-1));
    y = H*ρ_cur + randn(Nr);
    ΓIinv = diagm(ones(Cdouble,Nr));
    D = diagm(Nr-2,Nr,1 => 2ones(Cdouble,Nr-2), 0 => -ones(Cdouble,Nr-2) ,2 => -ones(Cdouble,Nr-2));
    yd = zeros(Cdouble,Nr-2);
    Γdinv = diagm(ones(Cdouble,Nr-2));

    # test the propsed sample 
    ρ_new,r_cp  = acceptSample(ρ_cur,ρ_prop,y,yd,ΓIinv,Γdinv,H,D);

    # results 
    ((ρ_new==ρ_cur) | (ρ_new==ρ_prop)) & (!isnan(r_cp)) & (!isinf(r_cp))
end


#TODO
function test_samplePosterior()
    true
end

function test_samplePosteriorMeanAndCov()
    true
end

function test_samplePosteriorModelMargin()
    true
end

function test_acceptSampleModelMargin()
    true
end


function test_samplePosteriorMargin()
    true
end

function test_acceptSampleMargin()
    true
end


@testset "Sampling" begin
    @test test_acceptSample()
    @test test_samplePosterior()
    @test test_samplePosteriorMeanAndCov()

    @test test_samplePosteriorModelMargin()
    @test test_acceptSampleModelMargin()

    @test test_samplePosteriorMargin()
    @test test_acceptSampleMargin()
end