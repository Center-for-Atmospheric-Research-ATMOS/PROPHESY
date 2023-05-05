
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


# Metropolis Hasting quadratic: MHquad.jl
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


function test_samplePosterior()
    # number of sample
    Ns = 5000;
    # dimension of each sample
    Nr = 10;
    ρ_GT  = 10ones(Cdouble,Nr);
    # measurement model
    H = diagm(0 => ones(Cdouble,Nr) ,1 => 0.5ones(Cdouble,Nr-1));
    # data
    y = H*ρ_GT + randn(Nr);
    # inverse covariance of the data noise
    ΓIinv = diagm(ones(Cdouble,Nr));
    # regularization operator
    D = diagm(Nr-2,Nr,1 => 2ones(Cdouble,Nr-2), 0 => -ones(Cdouble,Nr-2) ,2 => -ones(Cdouble,Nr-2));
    # a priori data
    yd = zeros(Cdouble,Nr-2);
    # regularization inverse covariance (strength and correlation of the regularization)
    Γdinv = diagm(ones(Cdouble,Nr-2));
    # sampling covariance matrix (used for generating samples in the communication mechanism)
    w = 0.25*0.1*ones(Cdouble,Nr); 
    Γsqrt = real(sqrt(corrCovariance(w;cor_len=5.0)));

    # draw the samples
    ρ_all,deltaU = samplePosterior(0.0ρ_GT,Γsqrt,y,yd,ΓIinv,Γdinv,H,D;Ns=Ns);
    
    # conditions
    cond1 = (((Ns+1,Nr)==size(ρ_all)) & (!isnan(sum(ρ_all))) & (!isinf(sum(ρ_all))))
    cond2 = ((Ns==length(deltaU)) & (!isnan(sum(deltaU))) & (!isinf(sum(deltaU))))

    # results
    cond1 & cond2
end

function test_samplePosteriorMeanAndCov()
    # number of sample and burn in period
    Ns    = 500000;
    Nburn = 1000;
    # dimension of each sample
    Nr = 10;
    ρ_GT  = 10ones(Cdouble,Nr);
    # measurement model
    H = diagm(0 => ones(Cdouble,Nr) ,1 => 0.5ones(Cdouble,Nr-1));
    # data
    ΓI = 0.1diagm(ones(Cdouble,Nr));
    y = H*ρ_GT + ΓI*randn(Nr);
    # inverse covariance of the data noise
    ΓIinv = inv(ΓI);
    # regularization operator
    D = diagm(Nr-2,Nr,1 => 2ones(Cdouble,Nr-2), 0 => -ones(Cdouble,Nr-2) ,2 => -ones(Cdouble,Nr-2));
    # a priori data
    yd = zeros(Cdouble,Nr-2);
    # regularization inverse covariance (strength and correlation of the regularization)
    Γdinv = diagm(ones(Cdouble,Nr-2));
    # sampling covariance matrix (used for generating samples in the communication mechanism)
    w = 0.25*0.1*ones(Cdouble,Nr); 
    Γsqrt = real(sqrt(corrCovariance(w;cor_len=5.0)));

    # mean and covariance estimate from samples
    μρ_cum,Γρ_cum,deltaU = samplePosteriorMeanAndCov(0.0ρ_GT,Γsqrt,y,yd,ΓIinv,Γdinv,H,D;Ns=Ns,Nburn=Nburn);

    # conditions
    cond1 = ((Nr==length(μρ_cum)) & (!isnan(sum(μρ_cum))) & (!isinf(sum(μρ_cum))))
    cond2 = (((Nr,Nr)==size(Γρ_cum)) & (!isnan(sum(Γρ_cum))) & (!isinf(sum(Γρ_cum))))
    cond3 = ((Ns==length(deltaU)) & (!isnan(sum(deltaU))) & (!isinf(sum(deltaU))))
    cond4 = ((sqrt(0.5*(ρ_GT-μρ_cum)'*inv(Γρ_cum)*(ρ_GT-μρ_cum))/Nr)<1.5) # 3σ: this test should almost never fail in practice for big enough sample size
    
    # results
    cond1 & cond2 & cond3 & cond4
end


function test_acceptSampleMargin()
    Nr = 10;
    ρ_cur  = 10ones(Cdouble,Nr);
    ρ_prop = 9ones(Cdouble,Nr);
    H = diagm(0 => ones(Cdouble,Nr) ,1 => 0.5ones(Cdouble,Nr-1));
    y = H*ρ_cur + randn(Nr);
    ΓIinv = diagm(ones(Cdouble,Nr));
    D = diagm(Nr-2,Nr,1 => 2ones(Cdouble,Nr-2), 0 => -ones(Cdouble,Nr-2) ,2 => -ones(Cdouble,Nr-2));
    yd = zeros(Cdouble,Nr-2);
    Γdinv = diagm(ones(Cdouble,Nr-2));
    ΓHΓyinv = 0.01diagm(ones(Cdouble,Nr));

    # check out the proposed sample
    ρ_new,r_cp = acceptSampleMargin(ρ_cur,ρ_prop,y,yd,ΓIinv,Γdinv,H,D,ΓHΓyinv);

    # results 
    ((ρ_new==ρ_cur) | (ρ_new==ρ_prop)) & (!isnan(r_cp)) & (!isinf(r_cp))
end

function test_samplePosteriorMargin()
    # number of sample
    Ns = 5000;
    # dimension of each sample
    Nr = 10;
    ρ_GT  = 10ones(Cdouble,Nr);
    # measurement model
    H = diagm(0 => ones(Cdouble,Nr) ,1 => 0.5ones(Cdouble,Nr-1));
    # data
    y = H*ρ_GT + randn(Nr);
    # inverse covariance of the data noise
    ΓIinv = diagm(ones(Cdouble,Nr));
    # normalized covariance matrix of the error in the measurement model
    ΓHΓyinv = 0.01diagm(ones(Cdouble,Nr));
    # regularization operator
    D = diagm(Nr-2,Nr,1 => 2ones(Cdouble,Nr-2), 0 => -ones(Cdouble,Nr-2) ,2 => -ones(Cdouble,Nr-2));
    # a priori data
    yd = zeros(Cdouble,Nr-2);
    # regularization inverse covariance (strength and correlation of the regularization)
    Γdinv = diagm(ones(Cdouble,Nr-2));
    # sampling covariance matrix (used for generating samples in the communication mechanism)
    w = 0.25*0.1*ones(Cdouble,Nr); 
    Γsqrt = real(sqrt(corrCovariance(w;cor_len=5.0)));

    # draw samples
    ρ_all,deltaU = samplePosteriorMargin(ρ_GT,Γsqrt,y,yd,ΓIinv,Γdinv,H,D,ΓHΓyinv;Ns=Ns);

    # conditions
    cond1 = (((Ns+1,Nr)==size(ρ_all)) & (!isnan(sum(ρ_all))) & (!isinf(sum(ρ_all))))
    cond2 = ((Ns==length(deltaU)) & (!isnan(sum(deltaU))) & (!isinf(sum(deltaU))))

    # results
    cond1 & cond2
end



# Metropolis Hasting quadratic with boundary conditions: MHquadBoundary.jl

function test_acceptSampleModelMargin()
    Nr = 10;
    ρ_cur  = 10ones(Cdouble,Nr);
    ρ_prop = 9ones(Cdouble,Nr);
    H = diagm(0 => ones(Cdouble,Nr) ,1 => 0.5ones(Cdouble,Nr-1));
    ΓH = zeros(Cdouble,Nr,Nr,Nr);
    [ΓH[i,i,i] = 0.01 for i in 1:Nr]
    y = H*ρ_cur + randn(Nr);
    ΓIinv = diagm(ones(Cdouble,Nr));
    D = diagm(Nr-2,Nr,1 => 2ones(Cdouble,Nr-2), 0 => -ones(Cdouble,Nr-2) ,2 => -ones(Cdouble,Nr-2));
    DΓinvD = D'diagm(ones(Cdouble,Nr-2))*D;
    ρB = [10.0; 10.0];
    σB = [0.1; 0.1];

    # check out proposed sample
    ρ_new,r_cp = acceptSampleModelMargin(ρ_cur,ρ_prop,y,ΓIinv,H,ΓH,DΓinvD,ρB,σB)

    # results 
    ((ρ_new==ρ_cur) | (ρ_new==ρ_prop)) & (!isnan(r_cp)) & (!isinf(r_cp))
end

function test_samplePosteriorModelMargin()
    # number of sample
    Ns = 5000;
    # dimension of each sample
    Nr = 10;
    ρ_GT  = 10ones(Cdouble,Nr);
    # measurement model
    H = diagm(0 => ones(Cdouble,Nr) ,1 => 0.5ones(Cdouble,Nr-1));
    # data
    y = H*ρ_GT + randn(Nr);
    # inverse covariance of the data noise
    ΓIinv = diagm(ones(Cdouble,Nr));
    # covariance matrixes of the error in the measurement model
    ΓH = zeros(Cdouble,Nr,Nr,Nr);
    [ΓH[i,i,i] = 0.01 for i in 1:Nr]
    # regularization operator
    D = diagm(Nr-2,Nr,1 => 2ones(Cdouble,Nr-2), 0 => -ones(Cdouble,Nr-2) ,2 => -ones(Cdouble,Nr-2));
    DΓinvD = D'diagm(ones(Cdouble,Nr-2))*D;
    # boundary values
    ρB = [10.0; 10.0];
    σB = [0.1; 0.1];
    # sampling covariance matrix (used for generating samples in the communication mechanism)
    w = 0.25*0.1*ones(Cdouble,Nr); 
    Γsqrt = real(sqrt(corrCovariance(w;cor_len=5.0)));

    # generate the chain of samples
    ρ_all = samplePosteriorModelMargin(ρ_GT,Γsqrt,y,ΓIinv,H,ΓH,DΓinvD,ρB,σB;Ns=Ns,psmooth=0.99);

    #results
    (((Ns+1,Nr)==size(ρ_all)) & (!isnan(sum(ρ_all))) & (!isinf(sum(ρ_all))))
end


function test_acceptSampleBoundary()
    Nr = 10;
    ρ_cur  = 10ones(Cdouble,Nr);
    ρ_prop = 9ones(Cdouble,Nr);
    H = diagm(0 => ones(Cdouble,Nr) ,1 => 0.5ones(Cdouble,Nr-1));
    ΓH = zeros(Cdouble,Nr,Nr,Nr);
    [ΓH[i,i,i] = 0.01 for i in 1:Nr]
    y = H*ρ_cur + randn(Nr);
    ΓIinv = diagm(ones(Cdouble,Nr));
    D = diagm(Nr-2,Nr,1 => 2ones(Cdouble,Nr-2), 0 => -ones(Cdouble,Nr-2) ,2 => -ones(Cdouble,Nr-2));
    DΓinvD = D'diagm(ones(Cdouble,Nr-2))*D;
    ρB = [10.0; 10.0];
    σB = [0.1; 0.1];

    # check out proposed sample
    ρ_new,r_cp = acceptSampleBoundary(ρ_cur,ρ_prop,y,ΓIinv,H,DΓinvD,ρB,σB);
    
    # results 
    ((ρ_new==ρ_cur) | (ρ_new==ρ_prop)) & (!isnan(r_cp)) & (!isinf(r_cp))
end

function test_samplePosteriorBoundary()
    # number of sample
    Ns = 5000;
    # dimension of each sample
    Nr = 10;
    ρ_GT  = 10ones(Cdouble,Nr);
    # measurement model
    H = diagm(0 => ones(Cdouble,Nr) ,1 => 0.5ones(Cdouble,Nr-1));
    # data
    y = H*ρ_GT + randn(Nr);
    # inverse covariance of the data noise
    ΓIinv = diagm(ones(Cdouble,Nr));
    # covariance matrixes of the error in the measurement model
    ΓH = zeros(Cdouble,Nr,Nr,Nr);
    [ΓH[i,i,i] = 0.01 for i in 1:Nr]
    # regularization operator
    D = diagm(Nr-2,Nr,1 => 2ones(Cdouble,Nr-2), 0 => -ones(Cdouble,Nr-2) ,2 => -ones(Cdouble,Nr-2));
    DΓinvD = D'diagm(ones(Cdouble,Nr-2))*D;
    # boundary values
    ρB = [10.0; 10.0];
    σB = [0.1; 0.1];
    # sampling covariance matrix (used for generating samples in the communication mechanism)
    w = 0.25*0.1*ones(Cdouble,Nr); 
    Γsqrt = real(sqrt(corrCovariance(w;cor_len=5.0)));

    # generate the chain of sample
    ρ_all = samplePosteriorBoundary(ρ_GT,Γsqrt,y,ΓIinv,H,DΓinvD,ρB,σB;Ns=Ns,psmooth=0.99);

    #results
    (((Ns+1,Nr)==size(ρ_all)) & (!isnan(sum(ρ_all))) & (!isinf(sum(ρ_all))))
end





@testset "Sampling" begin
    @test test_acceptSample()
    @test test_samplePosterior()
    @test test_samplePosteriorMeanAndCov()

    @test test_samplePosteriorModelMargin()
    @test test_acceptSampleModelMargin()

    @test test_samplePosteriorMargin()
    @test test_acceptSampleMargin()

    @test test_acceptSampleBoundary()
    @test test_samplePosteriorBoundary()
end