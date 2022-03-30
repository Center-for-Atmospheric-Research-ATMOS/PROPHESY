"""

# XPSmeas_simple.jl --
#
# XPSmeas_simple.jl implements a simplified version of the discretization
# of the different XPS configurations:
#   - planar interface
#   - microjet: cylindrical geometry
#   - particle/droplet: spherical model (TODO)
#
#------------------------------------------------------------------------------
#
# This file is part of the XPSmeas module which is licensed under the MIT "Expat" License:
#
# Copyright (C) 2021-2022,  Matthew Ozon.
#
#------------------------------------------------------------------------------

"""


##
## sharp edge approximation: geometrical distance
##
"""
    d_plane_P(x::Cdouble,y::Cdouble,z::Array{Cdouble,1},x0::Cdouble,y0::Cdouble,z0::Cdouble)

    In the planar geometry, e.g. flat surface spinning wheel, this function returns the distance
    from any point M:(x,y,z) in the sample to the surface in direction to some point P:(x0,y0,z0).
    Note that the planar interface is z=0 and that the sample sits in the volume z⩽0. Also, the
    distances are computed for fixed x and y.
"""
function d_plane_P(x::Cdouble,y::Cdouble,z::Array{Cdouble,1},x0::Cdouble,y0::Cdouble,z0::Cdouble)
    C = -z.*sqrt.(1.0.+((x0.-x).^2. +(y0.-y).^2)./((z0.-z).^2))
    C[C.<=0.0] .= 0.0;
    C
end


"""

    d_cylinder_P(r::Array{Cdouble,1},θ::Array{Cdouble,1},y::Cdouble,x0::Cdouble,y0::Cdouble,z0::Cdouble,μ0::Cdouble)

    In the cylindrical geometry, e.g. microjet, this function computes the distance from any
    point M:(r,θ,y) in the sample to the surface in direction to some point P:(x0,y0,z0).
    Note: M is given in cylindrical coordinates and P in Cartesian coordinates.
    Note: the mapping is computed only for one slice at y.

    The geometry of the sample is so that Oy is the symmetry axis (along the liquid stream direction),
    Oz is along the photon beam and Ox is the remaining axis in the orthogonal basis.
    θ is taken between Oz and the project of M onto the plane xOz and r = √(x^2+z^2)
    μ0 is the radius of the cylinder
"""
function d_cylinder_P(r::Array{Cdouble,1},θ::Array{Cdouble,1},y::Cdouble,x0::Cdouble,y0::Cdouble,z0::Cdouble,μ0::Cdouble)
    R0 = sqrt(z0^2+x0^2);
    θ0 = atan(x0,z0);
    A = -(R0*r*cos.(θ'.-θ0) .- r.^2).*sqrt.(R0^2 .- 2R0*r*cos.(θ'.-θ0) .+ r.^2 .+ (y0.-y).^2)./(R0^2 .- 2R0*r*cos.(θ'.-θ0) .+ r.^2)
    B = (1.0 .+ (y0-y).^2 ./(R0^2 .- 2R0*r*cos.(θ'.-θ0) .+ r.^2)).*(μ0^2 .- (R0^2*r.^2*(sin.(θ'.-θ0)).^2)./(R0^2 .- 2R0*r*cos.(θ'.-θ0) .+ r.^2))
    A[r.>μ0,:] .= 0.0
    B[r.>μ0,:] .= 0.0
    C = A+sqrt.(B);
    C[C.<0.0] .= 0.0;
    C
end

"""
    d_sphere_P(r::Array{Cdouble,1},θ::Array{Cdouble,1},φ::Array{Cdouble,1},x0::Cdouble,y0::Cdouble,z0::Cdouble,μ0::Cdouble)

    NOT IMPLEMENTED YET!

    In the spherical geometry, e.g. droplet, this function computes the distance from any
    point M:(r,θ,φ) in the sample to the surface in direction to some point P:(x0,y0,z0).
    Note: M is given in spherical coordinates and P in Cartesian coordinates.

    The geometry of the sample is so that Oy is the symmetry axis (along the droplet stream direction),
    Oz is along the photon beam and Ox is the remaining axis in the orthogonal basis.
    θ is taken between Oz and the project of M onto the plane xOz, φ between M and Oy, and r = √(x^2+y^2+z^2)
    μ0 is the radius of the cylinder
"""
function d_sphere_P(r::Array{Cdouble,1},θ::Array{Cdouble,1},φ::Array{Cdouble,1},x0::Cdouble,y0::Cdouble,z0::Cdouble,μ0::Cdouble)
    # TODO
end


##
## discretized models
##
"""
    plane_gain_H(x::Array{Cdouble,1},y::Array{Cdouble,1},z::Array{Cdouble,1},x0::Cdouble,y0::Cdouble,z0::Cdouble,λe::Cdouble)

    Compute the volume integrales (exact for piecewise linear gain)

    H_{n,j,k} = ∭ e_n(z)e_j(x)e_k(y) e^{d_P(M)/λe} dzdxdy

    The arrays x, y and z are the discretization subdivisions
    P:(x0,y0,z0) is the point in Cartesian coordinates used for the computation of the distance d_P
    λe is the attenuation length in the Beer-Lambert model

    Note: does not include the profile of the beam light, but it can be add easily if known
"""
function plane_gain_H(x::Array{Cdouble,1},y::Array{Cdouble,1},z::Array{Cdouble,1},x0::Cdouble,y0::Cdouble,z0::Cdouble,λe::Cdouble)
    # compute the elementary integrals
    Azn = 0.5*[z[2]-z[1]; z[3:end]-z[1:end-2]; z[end]-z[end-1]]; # dz
    Axj = 0.5*[x[2]-x[1]; x[3:end]-x[1:end-2]; x[end]-x[end-1]]; # dx
    Ayk = 0.5*[y[2]-y[1]; y[3:end]-y[1:end-2]; y[end]-Y[end-1]]; # dy

    #compute the model
    H_zxy = zeros(Cdouble,length(z),length(x),length(y));

    for j in 1:length(x)
        for k in 1:length(y)
            H_zxy[:,j,k] = exp.(-d_plane_P(x[j],y[k],z,x0,y0,z0)/λe)
        end
    end

    H_z  = zeros(Cdouble,length(z));
    for n in 1:length(z)
        H_z[n]  = Azn[n]*Axj'*H_z_far[n,:,:]*Ayk
    end

    H_z,H_zxy,Azn,Axj,Ayk
end

"""
    cylinder_gain_H(r::Array{Cdouble,1},θ::Array{Cdouble,1},y::Array{Cdouble,1},x0::Cdouble,y0::Cdouble,z0::Cdouble,μ0::Cdouble,λe::Cdouble)

    Compute the volume integrales (exact for piecewise linear gain)

    H_{n,j,k} = ∭ e_n(r)e_j(θ)e_k(y) e^{d_P(M)/λe} rdrdθdy

    The arrays r, θ and y are the discretization subdivisions
    P:(x0,y0,z0) is the point in Cartesian coordinates used for the computation of the distance d_P
    μ0 is the radius of the cylinder
    λe is the attenuation length in the Beer-Lambert model

    Note: does not include the profile of the beam light, but it can be add easily if known
"""
function cylinder_gain_H(r::Array{Cdouble,1},θ::Array{Cdouble,1},y::Array{Cdouble,1},x0::Cdouble,y0::Cdouble,z0::Cdouble,μ0::Cdouble,λe::Cdouble)
    # compute the elementary integrals
    Arn = 0.5*r.*[r[2]-r[1]; r[3:end]-r[1:end-2]; r[end]-r[end-1]]; # rdr
    Aθj = 0.5*[θ[2]-θ[1]; θ[3:end]-θ[1:end-2]; θ[end]-θ[end-1]];    # dθ
    Ayk = 0.5*[y[2]-y[1]; y[3:end]-y[1:end-2]; y[end]-Y[end-1]];    # dy

    #compute the model
    H_rθy = zeros(Cdouble,length(r),length(θ),length(y));
    for k in 1:length(y)
        H_rθy[:,:,k] = exp.(-d_cylinder_P(r,θ,y[k],x0,y0,z0,μ0)/λe)
    end

    H_r  = zeros(Cdouble,length(r));
    for n in 1:length(r)
        H_r[n] = Arn[n]*Aθj'*H_rθy[n,:,:]*Ayk
    end
    H_r,H_rθy,Arn,Aθj,Ayk
end
