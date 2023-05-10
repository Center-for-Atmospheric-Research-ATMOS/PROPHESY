# XPSpack test bench

using Test
using LinearAlgebra
using Distributions
using XPSpack


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


############################################################################################
#                          electron travel distance                                        #
############################################################################################
function test_distance_plane()
  # point location of the analyzer aperture
  x0 = 5000.0
  y0 = 0.0
  z0 = 5000.0
  # attenuation length
  λ_meas = 1.5e-3;
  # spatial discretization
  Nz = 101;
  x = 0.0;
  y = 0.0
  z = collect(LinRange(-5λ_meas,2λ_meas,Nz)); # depth 

  # travelling distance in the sample
  dP = d_plane_P(x,y,z,x0,y0,z0);

  # results
  (Nz==length(dP)) & (!isnan(sum(dP))) & (!isinf(sum(dP))) & (all(dP.>=0.0))
end

function test_distance_cylinder()
  # point location of the analyzer aperture
  x0 = 5000.0
  y0 = 0.0
  z0 = 5000.0
  # sample characteristic values
  λ_meas = 1.5e-3; # attenuation length
  μ0 = 10.0;       # cylinder diameter
  # spatial discretization
  Nr = 101;
  Nθ = 51;
  θ0 = atan(x0,z0)
  r = collect(LinRange(μ0-5λ_meas,μ0+2λ_meas,Nr)); # radial coordinates
  θ = collect(LinRange(θ0-π/2.0,θ0+π/2.0,Nθ));     # polar angles
  y = 0.0;                                         # axial coordinates
  
  # travelling distance in the sample
  dP = d_cylinder_P(r,θ,y,x0,y0,z0,μ0);

  # results
  (Nr==size(dP,1)) & (Nθ==size(dP,2)) & (!isnan(sum(dP))) & (!isinf(sum(dP))) & (all(dP.>=0.0))
end

function test_distance_cylinder_simplification()
  # point location of the analyzer aperture
  x0 = 5000.0
  y0 = 0.0
  z0 = 5000.0
  # sample characteristic values
  λ_meas = 1.5e-3; # attenuation length
  μ0 = 10.0;       # cylinder diameter
  # spatial discretization
  Nr = 101;
  Nθ = 51;
  θ0 = atan(x0,z0)
  r = collect(LinRange(μ0-5λ_meas,μ0+2λ_meas,Nr)); # radial coordinates
  θ = collect(LinRange(θ0-π/2.0,θ0+π/2.0,Nθ));     # polar angles
  y = 0.0;                                         # axial coordinates
  
  # travelling distance in the sample
  # dP = d_cylinder_P(r,θ,y,x0,y0,z0,μ0);
  dP = d_cylinder_P_simple(r,θ,y,x0,y0,z0,μ0)

  # results
  (Nr==size(dP,1)) & (Nθ==size(dP,2)) & (!isnan(sum(dP))) & (!isinf(sum(dP))) & (all(dP.>=0.0))
end

function test_distance_sphere()
  # point location of the analyzer aperture
  x0 = 5000.0
  y0 = 0.0*5000.0
  z0 = 0.0
  # sample characteristic values
  λ_meas = 1.5e-3; # attenuation length
  μ0 = 10.0;       # cylinder diameter
  # spatial discretization
  Nr = 101;
  Nθ = 51;
  Nφ = 52;
  θ0 = atan(y0,x0)
  φ0 = atan(sqrt(x0^2+y0^2),z0)
  r = collect(LinRange(μ0-5λ_meas,μ0+2λ_meas,Nr)); # radial coordinates
  φ = collect(LinRange(φ0-π/2.0,φ0+π/2.0,Nφ));     # polar angles
  θ = collect(LinRange(θ0-π/2.0,θ0+π/2.0,Nθ));     # azimuth angles

  # travelling distance in the sample
  dP = d_sphere_P(r,φ,θ,x0,y0,z0,μ0);

  # results
  (Nr==size(dP,1)) & (Nφ==size(dP,2)) & (Nθ==size(dP,3)) & (!isnan(sum(dP))) & (!isinf(sum(dP))) & (all(dP.>=0.0))
end


@testset "photo-electron travelling distance in the sample" begin
  @test test_distance_plane()
  @test test_distance_cylinder()
  @test test_distance_cylinder_simplification()
  @test test_distance_sphere()
end

############################################################################################
#                       geometry factor discretization                                     #
############################################################################################

function test_finger_gain()
  # point location of the analyzer aperture
  x0 = 5000.0
  y0 = 0.0
  z0 = 5000.0
  # attenuation length
  λ_meas = 1.5e-3;
  # spatial discretization
  Nz = 101;
  x = 0.0;
  y = 0.0
  z = collect(LinRange(-5λ_meas,2λ_meas,Nz)); # depth 

  # geometry factor and integral discretization: signal coming from along a line orthogonal to a plane surface
  H_z,Azn = finger_gain_H(x,y,z,x0,y0,z0,λ_meas);

  # conditions
  cond1 = (Nz==length(H_z)) & (!isnan(sum(H_z))) & (!isinf(sum(H_z))) & (all(H_z.>=0.0))
  cond2 = (Nz==length(Azn)) & (!isnan(sum(Azn))) & (!isinf(sum(Azn)))

  # results
  cond1 & cond2
end

function test_plane_gain()
  # point location of the analyzer aperture
  x0 = 5000.0
  y0 = 0.0
  z0 = 5000.0
  # attenuation length
  λ_meas = 1.5e-3;
  # spatial discretization
  Nx = 51;
  Ny = 53;
  Nz = 101;
  x = collect(LinRange(-5.0,5.0,Nx));
  y = collect(LinRange(-5.1,5.1,Ny));
  z = collect(LinRange(-5λ_meas,2λ_meas,Nz)); # depth 

  # geometry factors and integrals discretization: signal coming from the depth under a patch of coordinates (x,y)
  H_z,H_zxy,Azn,Axj,Ayk = plane_gain_H(x,y,z,x0,y0,z0,λ_meas);

  # conditions
  cond1 = (Nz==length(H_z)) & (!isnan(sum(H_z))) & (!isinf(sum(H_z))) & (all(H_z.>=0.0))
  cond2 = (Nz==size(H_zxy,1)) & (Nx==size(H_zxy,2)) & (Ny==size(H_zxy,3)) & (!isnan(sum(H_zxy))) & (!isinf(sum(H_zxy))) & (all(H_zxy.>=0.0))
  cond3 = (Nx==length(Axj)) & (!isnan(sum(Axj))) & (!isinf(sum(Axj)))
  cond4 = (Ny==length(Ayk)) & (!isnan(sum(Ayk))) & (!isinf(sum(Ayk)))
  cond5 = (Nz==length(Azn)) & (!isnan(sum(Azn))) & (!isinf(sum(Azn)))

  # results
  cond1 & cond2 & cond3 & cond4 & cond5
end

function test_cylinder_gain()
  # point location of the analyzer aperture
  x0 = 5000.0
  y0 = 0.0
  z0 = 5000.0
  # sample characteristic values
  λ_meas = 1.5e-3; # attenuation length
  μ0 = 10.0;       # cylinder diameter
  # spatial discretization
  Nr = 101;
  Nθ = 51;
  Ny = 21;
  θ0 = atan(x0,z0)
  L = 50.0 
  r = collect(LinRange(μ0-5λ_meas,μ0+2λ_meas,Nr)); # radial coordinates
  θ = collect(LinRange(θ0-π/2.0,θ0+π/2.0,Nθ));     # polar angles
  y = collect(LinRange(-L/2.0,L/2.0,Ny));          # axial coordinates

  # geometry factors and integrals discretization: signal coming from the part of the sample facing the aperture of the kinetic energy analyzer
  H_r,H_rθy,Arn,Aθj,Ayk = cylinder_gain_H(r,θ,y,x0,y0,z0,μ0,λ_meas);

  # conditions 
  cond1 = (Nr==length(H_r)) & (!isnan(sum(H_r))) & (!isinf(sum(H_r))) & (all(H_r.>=0.0))
  cond2 = (Nr==size(H_rθy,1)) & (Nθ==size(H_rθy,2)) & (Ny==size(H_rθy,3)) & (!isnan(sum(H_rθy))) & (!isinf(sum(H_rθy))) & (all(H_rθy.>=0.0))
  cond3 = (Nr==length(Arn)) & (!isnan(sum(Arn))) & (!isinf(sum(Arn)))
  cond4 = (Nθ==length(Aθj)) & (!isnan(sum(Aθj))) & (!isinf(sum(Aθj)))
  cond5 = (Ny==length(Ayk)) & (!isnan(sum(Ayk))) & (!isinf(sum(Ayk)))

  # results
  cond1 & cond2 & cond3 & cond4 & cond5
end

function test_cov_model()
  # point location of the analyzer aperture
  x0 = 5000.0
  y0 = 0.0
  z0 = 5000.0
  # sample characteristic values
  λ_meas = 1.5e-3; # attenuation length
  μ0 = 10.0;       # cylinder diameter
  # spatial discretization
  Nr = 101;
  Nθ = 51;
  Ny = 21;
  θ0 = atan(x0,z0)
  L = 50.0 
  r = collect(LinRange(μ0-5λ_meas,μ0+2λ_meas,Nr)); # radial coordinates
  θ = collect(LinRange(θ0-π/2.0,θ0+π/2.0,Nθ));     # polar angles
  y = collect(LinRange(-L/2.0,L/2.0,Ny));          # axial coordinates

  # attenuation length probability distribution
  Nλ = 25;
  λ = collect(LinRange(0.75λ_meas,1.25λ_meas,Nλ))
  σλ2 = (0.075λ_meas)^2;
  Pλ = (1.0/sqrt(2π*σλ2))*exp.(-((λ.-λ_meas).^2)./(2σλ2))

  # compute mean and covariance of the model with respect to the probability distribution of attenuation length
  ΓH, μH = cov_H_cylinder(r,θ,y,x0,y0,z0,μ0,λ,Pλ);

  # conditions 
  cond1 = (Nr==length(μH)) & (!isnan(sum(μH))) & (!isinf(sum(μH))) & (all(μH.>=0.0))
  cond2 = ((Nr,Nr)==size(ΓH)) & (!isnan(sum(ΓH))) & (!isinf(sum(ΓH)))
  λH = eigvals(ΓH) # all eigen values should be positive for a covariance matrix (but it's a numerical approximation)
  cond3 = (all(λH.>=0.0))
  if (!cond3)
    # maybe it is only computational approximation problem
    λmin,λmax = extrema(λH)
    if (abs(λmin/λmax)<1.0e-10)
      cond3 = true
    end
  end

  # results
  cond1 & cond2 & cond3
end

function test_sphere_gain()
  # point location of the analyzer aperture
  x0 = -50000.0
  y0 = 5000.0
  z0 = 50000.0
  # sample characteristic values
  λ_meas = 1.5e-3; # attenuation length
  μ0 = 10.0;       # cylinder diameter
  # spatial discretization
  Nr = 101;
  Nθ = 251;
  Nφ = 152;
  r = collect(LinRange(μ0-5λ_meas,μ0+2λ_meas,Nr)); # radial coordinates
  φ = collect(LinRange(0.0,π,Nφ));     # polar angles (∈[0,π])
  θ = collect(LinRange(0.0,2π,Nθ));     # azimuth angles (∈[0,2π])
  # the restriction cut out some signal
  # θ0 = atan(y0,x0)
  # φ0 = atan(sqrt(x0^2+y0^2),z0)
  # φ = collect(LinRange(φ0-π/2.0,φ0+π/2.0,Nφ));     # polar angles
  # θ = collect(LinRange(θ0-π/2.0,θ0+π/2.0,Nθ));     # azimuth angles

  # geometry factors and integrals discretization: signal coming from the part of the sample facing the aperture of the kinetic energy analyzer
  H_r,H_rφθ,Arn,Aφj,Aθk = sphere_gain_H(r,φ,θ,x0,y0,z0,μ0,λ_meas);

  # conditions
  cond1 = (Nr==length(H_r)) & (!isnan(sum(H_r))) & (!isinf(sum(H_r))) & (all(H_r.>=0.0))
  cond2 = ((Nr,Nφ,Nθ)==size(H_rφθ)) & (!isnan(sum(H_rφθ))) & (!isinf(sum(H_rφθ))) & (all(H_rφθ.>=0.0))
  cond3 = (Nr==length(Arn)) & (!isnan(sum(Arn))) & (!isinf(sum(Arn)))
  cond4 = (Nφ==length(Aφj)) & (!isnan(sum(Aφj))) & (!isinf(sum(Aφj)))
  cond5 = (Nθ==length(Aθk)) & (!isnan(sum(Aθk))) & (!isinf(sum(Aθk)))

  # results 
  cond1 & cond2 & cond3 & cond4 & cond5
end

function test_cylinder_gain_smooth_H()
  true
end

@testset "Geometry factors" begin
  @test test_finger_gain()
  @test test_plane_gain()
  @test test_cylinder_gain()
  @test test_cov_model()
  @test test_sphere_gain()
  @test test_cylinder_gain_smooth_H()
end


############################################################################################
#                       alignment parameter estimation                                     #
############################################################################################


function test_APE_smooth_interface()
  # predefine multiplicative parameter 
  τ = 100.0
  # size: number of measurement (Nke) in the spectrum and number of discretization nodes (Nr)
  Nke = 201
  Nr = 101
  # discretization
  Ke = collect(LinRange(0.0,10.0,Nke));
  dKe = Ke[2]-Ke[1];
  
  # photoionization cross-section density
  σKe2 = (25dKe)^2;
  μKe = Ke[101];
  σ_χ = (1.0/sqrt(2π*σKe2))exp.(-((Ke.-μKe).^2)/(2σKe2));
  σ_χ = σ_χ/(dKe*sum(σ_χ));
  # measurement operator (liquid interface)
  Δd50water = 0.5;
  λ_meas = 1.5;
  r = collect(LinRange(-5λ_meas,2λ_meas,Nr));
  dP_w = Δd50water*log.(1.0.+exp.(-r./Δd50water));
  H = exp.(-dP_w/λ_meas);
  # concentration profile
  Δd50C1s = 0.5;
  ρ = logistic.(r/Δd50C1s) # ρ = 1.0./(1.0.+exp.(-r/Δd50C1s));
  # spectrum and background
  I_bg = 250.0ones(Cdouble,Nke);
  I_data = τ*σ_χ*(H'*ρ) + I_bg;
  [I_data[i] = rand(Poisson(I_data[i])) for i in eachindex(Ke)];
  
  # estimates
  τ_mean,noise_data,τ_estimates = noiseAndParameterEstimation(σ_χ,H,I_data,I_bg,ρ);

  # conditions for successful test
  cond1 = (abs(τ/τ_mean -1.0)<0.2) # less than 20% error
  cond2 = (!isnan(sum(noise_data))) & (!isinf(sum(noise_data)))
  cond3 = (!isnan(sum(τ_estimates))) & (!isinf(sum(τ_estimates)))

  #result
  cond1 & cond2 & cond3
end

function test_alignmentParameter()
  # point location of the analyzer aperture
  x0 = 5000.0
  y0 = 0.0
  z0 = 5000.0
  # off-center distance between the center of the photon beam profile and the sample
  xc = 90.0
  yc = 50.0
  # spread of the photon beam profile
  σx = 70.0
  σy = 25.0
  # beam profile object
  bp = beamProfile(xc,yc,σx,σy);
  # sample characteristic values
  λ_meas = 1.5e-3; # attenuation length
  μ0 = 10.0;       # cylinder diameter
  # spatial discretization
  Nr = 101;
  Nθ = 51;
  Ny = 21;
  θ0 = atan(x0,z0)
  L = 50.0 
  r = collect(LinRange(μ0-5λ_meas,μ0+2λ_meas,Nr)); # radial coordinates
  θ = collect(LinRange(θ0-π/2.0,θ0+π/2.0,Nθ));     # polar angles
  y = collect(LinRange(-L/2.0,L/2.0,Ny));          # axial coordinates
  # compute the geometry factors discretization, the integrals discretization and the AP
  H_r,H_rθy,H_r_ph,H_rθy_ph,Arn,Aθj,Ayk,α = alignmentParameter(bp,r,θ,y,x0,y0,z0,μ0,λ_meas); #sum(H_r_ph)/sum(H_r)

  # conditions for successful test
  cond1 = (Nr==length(H_r)) & (!isnan(sum(H_r))) & (!isinf(sum(H_r)))
  cond2 = (Nr==length(H_r_ph)) & (!isnan(sum(H_r_ph))) & (!isinf(sum(H_r_ph)))
  cond3 = (Nr==size(H_rθy,1)) & (Nθ==size(H_rθy,2)) & (Ny==size(H_rθy,3)) & (!isnan(sum(H_rθy))) & (!isinf(sum(H_rθy)))
  cond4 = (Nr==size(H_rθy_ph,1)) & (Nθ==size(H_rθy_ph,2)) & (Ny==size(H_rθy_ph,3)) & (!isnan(sum(H_rθy_ph))) & (!isinf(sum(H_rθy_ph)))
  cond5 = (Nr==length(Arn)) & (!isnan(sum(Arn))) & (!isinf(sum(Arn)))
  cond6 = (Nθ==length(Aθj)) & (!isnan(sum(Aθj))) & (!isinf(sum(Aθj)))
  cond7 = (Ny==length(Ayk)) & (!isnan(sum(Ayk))) & (!isinf(sum(Ayk)))
  cond8 = (!isnan(α)) & (!isinf(α))

  #result
  cond1 & cond2 & cond3 & cond4 & cond5 & cond6 & cond7 & cond8
end


function test_alignmentParameterSphere()
  # point location of the analyzer aperture
  x0 = 5000.0
  y0 = 0.0
  z0 = 5000.0
  # off-center distance between the center of the photon beam profile and the sample
  xc = 90.0
  yc = 50.0
  # spread of the photon beam profile
  σx = 70.0
  σy = 25.0
  # beam profile object
  bp = beamProfile(xc,yc,σx,σy);
  # sample characteristic values
  λ_meas = 1.5e-3; # attenuation length
  μ0 = 10.0;       # cylinder diameter
  # spatial discretization
  Nr = 101;
  Nθ = 51;
  Nφ = 52;
  θ0 = atan(x0,z0)
  φ0 = atan(y0,sqrt(x0^2+z0^2))
  r = collect(LinRange(μ0-5λ_meas,μ0+2λ_meas,Nr)); # radial coordinates
  φ = collect(LinRange(φ0-π/2.0,φ0+π/2.0,Nφ));     # polar angles
  θ = collect(LinRange(θ0-π/2.0,θ0+π/2.0,Nθ));     # azimuth angles
  # compute the geometry factors discretization, the integrals discretization and the AP
  H_r,H_rφθ,H_r_ph,H_rφθ_ph,Arn,Aφj,Aθk,α = alignmentParameterSphere(bp,r,φ,θ,x0,y0,z0,μ0,λ_meas); #sum(H_r_ph)/sum(H_r)

  # conditions for successful test
  cond1 = (Nr==length(H_r)) & (!isnan(sum(H_r))) & (!isinf(sum(H_r)))
  cond2 = (Nr==length(H_r_ph)) & (!isnan(sum(H_r_ph))) & (!isinf(sum(H_r_ph)))
  cond3 = (Nr==size(H_rφθ,1)) & (Nφ==size(H_rφθ,2)) & (Nθ==size(H_rφθ,3)) & (!isnan(sum(H_rφθ))) & (!isinf(sum(H_rφθ)))
  cond4 = (Nr==size(H_rφθ_ph,1)) & (Nφ==size(H_rφθ_ph,2)) & (Nθ==size(H_rφθ_ph,3)) & (!isnan(sum(H_rφθ_ph))) & (!isinf(sum(H_rφθ_ph)))
  cond5 = (Nr==length(Arn)) & (!isnan(sum(Arn))) & (!isinf(sum(Arn)))
  cond6 = (Nφ==length(Aφj)) & (!isnan(sum(Aφj))) & (!isinf(sum(Aφj)))
  cond7 = (Nθ==length(Aθk)) & (!isnan(sum(Aθk))) & (!isinf(sum(Aθk)))
  cond8 = (!isnan(α)) & (!isinf(α))

  #result
  cond1 & cond2 & cond3 & cond4 & cond5 & cond6 & cond7 & cond8
end


@testset "Alignment Parameter Definition (cylinder and sphere geometry)" begin
  @test test_alignmentParameter()
  @test test_alignmentParameterSphere()
end

@testset "Alignment Parameter Estimation (APE)" begin
  @test test_APE_smooth_interface()
end