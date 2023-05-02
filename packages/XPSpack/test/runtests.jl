# XPSpack test bench

using Test
using LinearAlgebra
using XPSpack


#=
remaining functions to be tested

logistic

# possible geometry of the sample
fingerGeom, planeGeom, cylinderGeom

# distance and geometry factors
d_plane_P, d_sphere_P, plane_gain_H, finger_gain_H
cylinder_gain_H,alignmentParameter, d_cylinder_P, d_cylinder_P_simple
cov_H_cylinder
sphere_gain_H,d_sphere_P,alignmentParameterSphere
beamProfile

# function for noise estimation
noiseAndParameterEstimation
=#

############################################################################################
# intrinsic properties of the sample: photoionization cross-section and attenuation length #
############################################################################################

function test_λe()
  all(λe.(collect(250.0:100.0:2000.0)).>=0.0)
end

function test_λe_exp()
  all(λe_exp.(collect(250.0:100.0:2000.0)).>=0.0)
end

function test_σ_C1s_exp()
  all(σ_C1s_exp.(collect(200.0:100.0:2000.0)).>=0.0)
end

function test_σ_O1s_exp()
  all(σ_O1s_exp.(collect(200.0:100.0:2000.0)).>=0.0)
end

function test_σ_S2p_exp()
  all(σ_S2p_exp.(collect(200.0:100.0:2000.0)).>=0.0)
end

function test_σ_cs_orb()
  ((σ_cs_orb(1000.0,"C1s")>=0.0) & (σ_cs_orb(1000.0,"O1s")>=0.0) & (σ_cs_orb(1000.0,"S2p")>=0.0))
end

function test_σ_bg()
  all(σ_bg.(collect(200.0:100.0:2000.0)).>=0.0)
end

function test_σ_bg_density()
  all(σ_bg_density(collect(200.0:100.0:2000.0),1000.0,5.0).>=0.0)
end



@testset "Photionization cross-sections and electron attenuation length" begin
  @test test_λe()
  @test test_λe_exp()
  @test test_σ_C1s_exp()
  @test test_σ_O1s_exp()
  @test test_σ_S2p_exp()
  @test test_σ_cs_orb()
  @test test_σ_bg()
  @test test_σ_bg_density()
end


############################################################################################
#                      peak fitting and background removal                                 #
############################################################################################
function test_baseline_removal()
  # dummy data
  Ny = 200
  y = zeros(Cdouble,Ny) 
  # regularization matrix and parameters
  D_2nd = diagm(Ny-2,Ny,1 => 2ones(Cdouble,Ny-2), 0 => -ones(Cdouble,Ny-2) ,2 => -ones(Cdouble,Ny-2));
  κb = 0.01;
  λb = 1.0e5;
  # estimate background
  z = baseline_removal(y,λb,κb,D_2nd);
  # result
  abs(sum(z))<=sqrt(Ny)
end

function test_EM_peaks()
  τm = [1.0/3.0; 1.0/3.0; 1.0/3.0];
  μm = [290.3; 291.9; 293.5];
  σm = [290.3; 291.9; 293.5]/500.0;
  Np = length(τm);
  τt = zeros(Cdouble,Np);
  μt = zeros(Cdouble,Np);
  σt = zeros(Cdouble,Np);

  Ke = collect(286.0:0.05:298.0);
  μ_XR = (τm[1]/sqrt(2π*σm[1]))*exp.(-((Ke.-μm[1]).^2)./(2.0σm[1]^2)) .+ (τm[2]/sqrt(2π*σm[2]))*exp.(-((Ke.-μm[2]).^2)./(2.0σm[2]^2)) .+ (τm[3]/sqrt(2π*σm[3]))*exp.(-((Ke.-μm[3]).^2)./(2.0σm[3]^2))

  # estimate the peaks centers and spreads
  Niter = 100
  τt,μt,σt = EM_peaks(Ke,μ_XR,τm,μm,σm,Niter)
  # result 
  ((norm(μt.-μm)<1.0) & (norm(σt.-σm)<1.0) & (norm(τt-τm)<1.0))
end

function test_cross_section_spread_function()
  # dummy data
  σnoise = 0.001;
  τm = [1.0/3.0; 1.0/3.0; 1.0/3.0];
  μm = [290.3; 291.9; 293.5];
  σm = [290.3; 291.9; 293.5]/500.0;
  Ke = collect(286.0:0.05:298.0);
  μ_XR = (τm[1]/sqrt(2π*σm[1]))*exp.(-((Ke.-μm[1]).^2)./(2.0σm[1]^2)) .+ (τm[2]/sqrt(2π*σm[2]))*exp.(-((Ke.-μm[2]).^2)./(2.0σm[2]^2)) .+ (τm[3]/sqrt(2π*σm[3]))*exp.(-((Ke.-μm[3]).^2)./(2.0σm[3]^2))

  # optional parameters
  σ_2nd=0.001;
  Nbfgs = 10;
  Nsearch = 10;

  # estimate the cross-section density
  Xend,Hend,_,Nlast,R,σ_R = cross_section_spread_function(μ_XR,Ke,σnoise;Nbfgs=Nbfgs,Nsearch=Nsearch,σ_2nd=σ_2nd)

  # results 
  (!isnan(R)) & (!isinf(R)) & (!isnan(σ_R)) & (!isinf(σ_R)) & (!isnan(Nlast)) & (!isinf(Nlast)) & (!isnan(sum(Xend))) & (!isinf(sum(Xend))) & (!isnan(sum(Hend))) & (!isinf(sum(Hend)))
end

function test_cross_section_spread_function_sample()
  # dummy data
  σnoise = 0.001;
  τm = [1.0/3.0; 1.0/3.0; 1.0/3.0];
  μm = [290.3; 291.9; 293.5];
  σm = [290.3; 291.9; 293.5]/500.0;
  Ke = collect(286.0:0.05:298.0);
  μ_XR = (τm[1]/sqrt(2π*σm[1]))*exp.(-((Ke.-μm[1]).^2)./(2.0σm[1]^2)) .+ (τm[2]/sqrt(2π*σm[2]))*exp.(-((Ke.-μm[2]).^2)./(2.0σm[2]^2)) .+ (τm[3]/sqrt(2π*σm[3]))*exp.(-((Ke.-μm[3]).^2)./(2.0σm[3]^2))

  # optional parameters
  σ_2nd=0.001;
  Nbfgs = 10;
  Nsearch = 10;
  Nloop  = 10;

  # estimate the cross-section density
  Xend1,μ_XR1,Γ_XR1,R1,σ_R1,R_samples1 = cross_section_spread_function_sample(μ_XR.+σnoise*randn(length(Ke)),Ke,σnoise;Nbfgs=Nbfgs,Nsearch=Nsearch,N_sample=Nloop,σ_2nd=σ_2nd);

  # results: make sure all the variables are not NaN or Inf (Xend1,μ_XR1,Γ_XR1,R1,σ_R1,R_samples1)
  (length(R_samples1)==(Nloop+1)) & (size(Xend1,1)==length(Ke)) & (size(Xend1,2)==(Nloop+1)) & (length(μ_XR1)==length(Ke)) & (!isnan(R1)) & (!isinf(R1)) & (!isnan(σ_R1)) & (!isinf(σ_R1)) & (!isnan(sum(R_samples1))) & (!isinf(sum(R_samples1))) & (!isnan(sum(Xend1))) & (!isinf(sum(Xend1))) & (!isnan(sum(μ_XR1))) & (!isinf(sum(μ_XR1))) & (!isnan(sum(Γ_XR1))) & (!isinf(sum(Γ_XR1)))
end



@testset "Cross-section density and peak fitting" begin
  @test test_baseline_removal()
  @test test_EM_peaks()
  @test test_cross_section_spread_function()
  @test test_cross_section_spread_function_sample()
end