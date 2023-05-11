
#------------------------------------------------------------------------------
#
# This file is part of the XPSpack module which is licensed under CC-BY 4.0 license.
#
# Copyright (C) 2022,  Matthew Ozon, Konstantin Tumashevich and Nønne L. Prisle.
#
#------------------------------------------------------------------------------

"""
    beamProfile: structure describing the extent of the spread in space of the photon beam
    The beam density is taken to be a cross section of the parallel beam (orthogonal to the propagation direction)
    (xc,yc): coordinate of the center
    (σx,σy): spatial spread in x and y directions
"""
mutable struct beamProfile
    # coordinates of the center of the profile
    xc::Cdouble
    yc::Cdouble

    # spread of the beam in x and y direction 
    σx::Cdouble
    σy::Cdouble

    function beamProfile()
        new(0.0,0.0,1.0,1.0)
    end
    function beamProfile(xc_::Cdouble,yc_::Cdouble,σx_::Cdouble,σy_::Cdouble)
        new(xc_,yc_,σx_,σy_)
    end
    function beamProfile(ws::beamProfile)
        new(ws.xc,ws.yc,ws.σx,ws.σy)
    end
end

"""
    alignmentParameter(bp::beamProfile,r::Array{Cdouble,1},θ::Array{Cdouble,1},y::Array{Cdouble,1},x0::Cdouble,y0::Cdouble,z0::Cdouble,μ0::Cdouble,λe::Cdouble)

    Compute the volume integrales (exact for piecewise linear gain)

    H_{n,j,k} = ∭ e_n(r)e_j(θ)e_k(y) f_beam(r,θ,y) e^{-d_P(M)/λe} rdrdθdy
    H_{n,j,k} = ∭ e_n(r)e_j(θ)e_k(y) e^{-d_P(M)/λe} rdrdθdy

    The arrays r, θ and y are the discretization subdivisions
    P:(x0,y0,z0) is the point in Cartesian coordinates used for the computation of the distance d_P
    μ0 is the radius of the cylinder
    λe is the attenuation length in the Beer-Lambert model
    bp: a beam profile object that describes the deviation from the center of the target and the spatial spread of the beam

    Note: does not include the profile of the beam light, 
    but it can be add easily if known, see alignmentParameter
"""
function alignmentParameter(bp::beamProfile,r::Array{Cdouble,1},θ::Array{Cdouble,1},y::Array{Cdouble,1},x0::Cdouble,y0::Cdouble,z0::Cdouble,μ0::Cdouble,λe::Cdouble)
    # compute the elementary integrals
    Arn = 0.5*r.*[r[2]-r[1]; r[3:end]-r[1:end-2]; r[end]-r[end-1]]; # rdr
    Aθj = 0.5*[θ[2]-θ[1]; θ[3:end]-θ[1:end-2]; θ[end]-θ[end-1]];    # dθ
    Ayk = 0.5*[y[2]-y[1]; y[3:end]-y[1:end-2]; y[end]-y[end-1]];    # dy

    #compute the model
    H_rθy    = zeros(Cdouble,length(r),length(θ),length(y));
    H_rθy_ph = zeros(Cdouble,length(r),length(θ),length(y));
    bpx      = (1.0/(2π*bp.σx*bp.σy))*exp.(-((r*sin.(θ)'.-bp.xc).^2)./(2.0*bp.σx^2));
    for k in 1:length(y)
        H_rθy[:,:,k] = exp.(-d_cylinder_P(r,θ,y[k],x0,y0,z0,μ0)/λe)
        f_p = bpx*exp.(-((y[k].-bp.yc).^2)./(2.0*bp.σy^2))'
        H_rθy_ph[:,:,k] = H_rθy[:,:,k].*f_p;
    end

    # radial model with and without beam profile
    H_r     = zeros(Cdouble,length(r));
    H_r_ph  = zeros(Cdouble,length(r));
    for n in 1:length(r)
        H_r[n]    = Arn[n]*Aθj'*H_rθy[n,:,:]*Ayk
        H_r_ph[n] = Arn[n]*Aθj'*H_rθy_ph[n,:,:]*Ayk
    end
    H_r,H_rθy,H_r_ph,H_rθy_ph,Arn,Aθj,Ayk,sum(H_r_ph)/sum(H_r)
end


"""
    alignmentParameterSphere(bp::beamProfile,r::Array{Cdouble,1},φ::Array{Cdouble,1},θ::Array{Cdouble,1},x0::Cdouble,y0::Cdouble,z0::Cdouble,μ0::Cdouble,λe::Cdouble)

    Compute the volume integrales (exact for piecewise linear gain)

    H_{n,j,k} = ∭ e_n(r)e_j(φ)e_k(θ) f_beam(r,θ,y) e^{-d_P(M)/λe} r^2sin(φ)drdφdθ
    H_{n,j,k} = ∭ e_n(r)e_j(φ)e_k(θ) e^{-d_P(M)/λe} r^2sin(φ)drdφdθ

    The geometry of the sample is so that Oy is the symmetry axis (along the droplet stream direction),
    Oz is along the photon beam and Ox is the remaining axis in the orthogonal basis.

    The arrays r, φ and θ are the discretization subdivisions (r: radial distance, φ polar angle, and θ azimuthal angle).
    The polar axis is Oy, φ is the polar angle between M and Oy, θ is taken between Oz and the project of M onto the plane xOz, and r = √(x^2+y^2+z^2)

    P:(x0,y0,z0) is the point in Cartesian coordinates used for the computation of the distance d_P
    μ0 is the radius of the cylinder
    λe is the attenuation length in the Beer-Lambert model
    bp: a beam profile object that describes the deviation from the center of the target and the spatial spread of the beam

    Note: does not include the profile of the beam light, 
    but it can be add easily if known, see alignmentParameterSphere
"""
function alignmentParameterSphere(bp::beamProfile,r::Array{Cdouble,1},φ::Array{Cdouble,1},θ::Array{Cdouble,1},x0::Cdouble,y0::Cdouble,z0::Cdouble,μ0::Cdouble,λe::Cdouble)
    # compute the elementary integrals
    Ar1 = (1.0/(r[2]-r[1]))*( (r[2]/3.0)*(r[2]^3-r[1]^3) - (1.0/4.0)*(r[2]^4-r[1]^4) );
    ArN = (1.0/(r[end]-r[end-1]))*( (1.0/4.0)*(r[end]^4-r[end-1]^4) - (r[end-1]/3.0)*(r[end]^3-r[end-1]^3));
    Arn = [Ar1; (1.0/12.0)*((r[3:end].^3-r[1:end-2].^3) .+ r[2:end-1].*(r[3:end].^2-r[1:end-2].^2) .+ (r[3:end]-r[1:end-2])); ArN]; # r^2dr

    Aφ1 = cos(φ[1]) - (1.0/(φ[2]-φ[1]))*(sin(φ[2]) - sin(φ[1]));
    AφJ = (1.0/(φ[end]-φ[end-1]))*(sin(φ[end]) - sin(φ[end-1])) - cos(φ[end]);
    Aφj = [Aφ1; ((sin.(φ[2:end-1])-sin.(φ[1:end-2]))./(φ[2:end-1]-φ[1:end-2]) - (sin.(φ[3:end])-sin.(φ[2:end-1]))./(φ[3:end]-φ[2:end-1])); AφJ];    # sinφdφ

    Aθk = 0.5*[θ[2]-θ[1]; θ[3:end]-θ[1:end-2]; θ[end]-θ[end-1]];    # dθ

    #compute the model
    H_rφθ    = zeros(Cdouble,length(r),length(φ),length(θ));
    H_rφθ_ph = zeros(Cdouble,length(r),length(φ),length(θ));
    y = r*cos.(φ');
    x = r*sin.(φ');
    
    xym = r*sin.(φ')
    for k in 1:length(θ)
        H_rφθ[:,:,k] = exp.(-d_sphere_P(r,φ,θ[k],x0,y0,z0,μ0)/λe)
        local f_p = (1.0/(2π*bp.σx*bp.σy)).*exp.(-((cos(θ[k])*xym.-bp.xc).^2/(2.0*bp.σx^2))).*exp.(-((sin(θ[k])*xym.-bp.yc).^2/(2.0*bp.σy^2)))
        H_rφθ_ph[:,:,k] = H_rφθ[:,:,k].*f_p;
    end

    # radial model with and without beam profile
    H_r     = zeros(Cdouble,length(r));
    H_r_ph  = zeros(Cdouble,length(r));
    for n in 1:length(r)
        H_r[n]    = Arn[n]*Aφj'*H_rφθ[n,:,:]*Aθk
        H_r_ph[n] = Arn[n]*Aφj'*H_rφθ_ph[n,:,:]*Aθk
    end
    H_r,H_rφθ,H_r_ph,H_rφθ_ph,Arn,Aφj,Aθk,sum(H_r_ph)/sum(H_r)
end


"""
   noiseAndParameterEstimation(σ_χ::Array{Cdouble,1},H::Array{Cdouble,1},I_data::Array{Cdouble,1},I_bg::Array{Cdouble,1},ρ::Array{Cdouble,1})

   returns the mean estimate of the multiplicative factor in the measurement model (including the alignment factor) and the samples τ_{ℓ,k} used to compute the mean
   and an estimate of the noise perturbations.

   σ_χ is the cross section density
   H is the geometrical factor
   I_data is the raw data (including the background)
   I_bg is the background
   ρ is a rough profile (potentially the bulk concentration in the volume and 0 outside)
"""
function noiseAndParameterEstimation(σ_χ::Array{Cdouble,1},H::Array{Cdouble,1},I_data::Array{Cdouble,1},I_bg::Array{Cdouble,1},ρ::Array{Cdouble,1})
   # estimate the noise using the measurement model (this could be computed from curve fits directly)
   F  = svd(σ_χ*H')
   noise_data = F.U[:,2:end]*(F.U[:,2:end]'*(I_data-I_bg));
   # estimate the signal of interest
   σ_data = I_data-(I_bg+noise_data)
   # compute the multiplicative factor estimates τ_{ℓ,k}
   τ_estimates = σ_data[σ_χ.>0.1*maximum(σ_χ)]./(σ_χ[σ_χ.>0.1*maximum(σ_χ)]*(H'*ρ));

   # return the mean value of the estimates, the estimated noise and the estimates
   mean(τ_estimates),noise_data,τ_estimates
end