#
# XPSmeas.jl --
#
# XPSmeas.jl is a module aiming at modelling XPS measurement, e.g. simulation
# of PE spectra, and inverting data, e.g. concentration profiles.
#
#
#------------------------------------------------------------------------------
#
# This file is part of the XPSpack module which is licensed under the MIT "Expat" License:
#
# Copyright (C) 2021-2022,  Matthew Ozon.
#
#------------------------------------------------------------------------------
"""

# XPSmeas.jl --
#
# XPSmeas.jl implements a simplified version of the discretization
# of the different XPS experiment configurations:
#   - planar interface
#   - microjet: cylindrical geometry
#   - particle/droplet: spherical model 
#
#------------------------------------------------------------------------------
#
# This file is part of the XPSpack module which is licensed under the MIT "Expat" License:
#
# Copyright (C) 2021-2022,  Matthew Ozon.
#
#------------------------------------------------------------------------------

"""


# the geometry structures assumes that
#  - the photon beam is along the z axis
#  - the liquid microjet is along the y axis (or the droplet axis motion, or, if no clear axis arise naturally, any direction orthogonal to the z axis)
#  - the x axis is the remaining axis that makes Oxyz an orthogonal frame of reference
#  - the reference of axis, namely O, is taken either at the center or on the surface of the sample
# maybe but not sure: the model derived from the geomStructs assumes that the incident light is orthogonal to the surface of the sample

# finger struct
"""
    fingerGeom is a mutable structure that contains the necessary parameters to describe a 1D geometry 
    where some signal may come from, hence the name finger.

    - x0,y0,z0: location of the analyzer's apperture
    - x,y:      coordinate of the electron source at the surface of the planar sample
    - z:        array of depth into the sample

    fingerGeom(x0_::Cdouble,y0_::Cdouble,z0_::Cdouble,x_::Cdouble,y_::Cdouble,z_::Array{Cdouble,1})
    creates an object of type fingerGeom
"""
mutable struct fingerGeom
    # coordinates of the analyzer's apperture (P)
    x0::Cdouble
    y0::Cdouble
    z0::Cdouble

    # coordinate of the finger at the sample's surface
    x::Cdouble
    y::Cdouble

    # probing depths along the finger (since the sampling can be non-uniform, I decide to keep an array instead of a range)
    z::Array{Cdouble,1}

    function fingerGeom()
        new(0.0,0.0,0.0, 0.0,0.0,Array{Cdouble,1}(undef,0))
    end
    function fingerGeom(x0_::Cdouble,y0_::Cdouble,z0_::Cdouble,x_::Cdouble,y_::Cdouble,z_::Array{Cdouble,1})
        new(x0_,y0_,z0_,x_,y_,z_)
    end
    function fingerGeom(ws::fingerGeom)
        new(ws.x0,ws.y0,ws.z0,ws.x,ws.y,ws.z)
    end
end

# plane struct: assumes a rectangular cuboid discretization volume
"""
    planeGeom is a mutable structure that contains the necessary parameters to describe a 3D geometry 
    where some signal may come from (some cuboid volume underneath the surface of a planar sample)

    - x0,y0,z0: location of the analyzer's apperture
    - x,y:      coordinates of the electron source at the surface of the planar sample
    - z:        array of depth into the sample

    planeGeom(x0_::Cdouble,y0_::Cdouble,z0_::Cdouble,x_::Array{Cdouble,1},y_::Array{Cdouble,1},z_::Array{Cdouble,1})
    creates an object of type planeGeom
"""
mutable struct planeGeom
    # coordinates of the analyzer's apperture (P)
    x0::Cdouble
    y0::Cdouble
    z0::Cdouble

    # coordinates of the illuminated area at the surface of the sample (the covered area may be bigger or smaller than the sample's dimension: the derived model should deal with that with the beam-light's profile)
    x::Array{Cdouble,1}
    y::Array{Cdouble,1}

    # probing depths
    z::Array{Cdouble,1}

    function planeGeom()
        new(0.0,0.0,0.0, Array{Cdouble,1}(undef,0),Array{Cdouble,1}(undef,0),Array{Cdouble,1}(undef,0))
    end
    function planeGeom(x0_::Cdouble,y0_::Cdouble,z0_::Cdouble,x_::Array{Cdouble,1},y_::Array{Cdouble,1},z_::Array{Cdouble,1})
        new(x0_,y0_,z0_,x_,y_,z_)
    end
    function planeGeom(ws::planeGeom)
        new(ws.x0,ws.y0,ws.z0,ws.x,ws.y,ws.z)
    end
end

# cylinder struct
"""
    cylinderGeom is a mutable structure that contains the necessary parameters to describe a 3D geometry 
    where some signal may come from (some cuboid volume underneath the surface of a planar sample)

    - x0,y0,z0: location of the analyzer's apperture (Cartesian coordinate)
    - μ0:       radius of the liquid microjet
    - θ,y:      coordinates of the electron source at the surface of the cylindrical sample (cylindrical coordinates)
    - r:        array of depth into the sample (radius in the sample)

    cylinderGeom(x0_::Cdouble,y0_::Cdouble,z0_::Cdouble,μ0_::Cdouble,r_::Array{Cdouble,1},θ_::Array{Cdouble,1},y_::Array{Cdouble,1})
    creates an object of type planeGeom
"""
mutable struct cylinderGeom
    # coordinates of the analyzer's apperture (P)
    x0::Cdouble
    y0::Cdouble
    z0::Cdouble

    # radius of the sample
    μ0::Cdouble 

    # probing depths (distance from the symmetry axis of the sample)
    r::Array{Cdouble,1}

    # coordinates of the illuminated area at the surface of the sample in polar coordinates (the covered area may be bigger or smaller than the sample's dimension: the derived model should deal with that with the beam-light's profile)
    θ::Array{Cdouble,1}
    y::Array{Cdouble,1}

    function cylinderGeom()
        new(0.0,0.0,0.0,0.0, Array{Cdouble,1}(undef,0),Array{Cdouble,1}(undef,0),Array{Cdouble,1}(undef,0))
    end
    function cylinderGeom(x0_::Cdouble,y0_::Cdouble,z0_::Cdouble,μ0_::Cdouble,r_::Array{Cdouble,1},θ_::Array{Cdouble,1},y_::Array{Cdouble,1})
        new(x0_,y0_,z0_,μ0_,r_,θ_,y_)
    end
    function cylinderGeom(ws::cylinderGeom)
        new(ws.x0,ws.y0,ws.z0,ws.μ0,ws.r,ws.θ,ws.y)
    end
end

# # acquisition struct: for each photon energy, this structure gather the information needed by the model. Create a dictionary of XPSacq that gathers all the photon energy to create the complete model
# """
#     XPSacq is a mutable structure that gathers in one place the information of an XPS acquisition (from a model point of view)

#     - ħν: photon energy
#     - μKe: analyzer's reference kinetic energy
#     - α: apparent solid angle (analyzer's apperture from the sample)
#     - T: transmission factor (depends mostly on μKe)
#     - Fν: photon flux
#     - Ke: array of probed kinetic energies
#     - σν: kinetic-energy-density-differential cross-section dσ/dKe → σ_{nl}(ħν) = ∫ dσ/dKe dKe
#     - λe: attenuation length (electron in the sample: Beer-Lambert's law)
#     - nl: orbitale's name 
# """
# mutable struct XPSacq
#     # main setup number
#     ħν::Cdouble                    # photon energy 
#     μKe::Cdouble                   # mean kinetic energies (where the anaylizer is measuring in the )

#     # analyzer
#     α::Cdouble                     # the apparent solid angle that the analyzer's apperture represents from the sample point of view
#     T::Cdouble                     # mystic gain of the analyzer TODO: demystify it (see robert's email)

#     # sampled kinetic energies
#     Fν::Cdouble                    # photon flux at the given photon energy
#     Ke::Array{Cdouble,1}           # kinetic-energy spectrum sampling points
#     σν::Array{Cdouble,1}           # kinetic-energy-density-differential cross-section dσ/dKe → σ_{nl}(ħν) = ∫ dσ/dKe dKe

#     # attenuation length
#     λe::Cdouble                    # 1D array for the electron attenuation length (w.r.t. K_e)

#     # geometry structure (WARNING: does it make sense to have it in this structure knowing that it's gonna be shared across many instance of the XPSacq, maybe leave it for the global structure XPSexp)
#     # tag for the orbitale (the differential cross section should already contain the orbitale information: energy-and-angula spectrum)
#     nl::String                     # orbitale

#     # default ctor (it is not really meaningful, it's more for the sake of having a default constructor)
#     function XPSacq() #
#         new(0.0,0.0,0.0,0.0,0.0,Array{Cdouble,1}(undef,0),Array{Cdouble,1}(undef,0),0.0,"C1s")
#     end

#     # ctor: this is physically relevant (as much as the meaning of the parameters used for the model)
#     function XPSacq(ħν_::Cdouble,μKe_::Cdouble,α_::Cdouble,T_::Cdouble,Fν_::Cdouble,Ke_::Array{Cdouble,1},σν_::Array{Cdouble,1},λe_::Cdouble;nl_::String="C1s")
#         if (length(Ke_)!=length(σν_))
#             throw("XPSacq: kinetic energy array and cross section array are not of the same length")
#         end
#         new(ħν_,μKe_,α_,T_,Fν_,Ke_,σν_,λe_,nl_)
#     end

#     # cptor
#     function XPSacq(ws::XPSacq) #
#         new(ws.ħν,ws.μKe,ws.α,ws.T,ws.Fν,ws.Ke,ws.σν,ws.λe,ws.nl)
#     end
# end


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
    d_cylinder_P_simple(r::Array{Cdouble,1},θ::Array{Cdouble,1},y::Cdouble,x0::Cdouble,y0::Cdouble,z0::Cdouble,μ0::Cdouble)

    Simplification of d_cylinder_P in the case μ0≪√((x0-x)^2+(z0-z)^2) and  |y0-y|≪√((x0-x)^2+(z0-z)^2)
"""
function d_cylinder_P_simple(r::Array{Cdouble,1},θ::Array{Cdouble,1},y::Cdouble,x0::Cdouble,y0::Cdouble,z0::Cdouble,μ0::Cdouble)
    θ0 = atan(x0,z0)
    A = -r*cos.(θ'.-θ0)
    B = (μ0^2 .- (r.^2*(sin.(θ'.-θ0)).^2))
    A[r.>μ0,:] .= 0.0
    B[r.>μ0,:] .= 0.0
    C = A+sqrt.(B);
    C[C.<0.0] .= 0.0;
    C
end

"""
    d_sphere_P(r::Array{Cdouble,1},θ::Array{Cdouble,1},φ::Array{Cdouble,1},x0::Cdouble,y0::Cdouble,z0::Cdouble,μ0::Cdouble)

    In the spherical geometry, e.g. droplet, this function computes the distance from any
    point M:(r,θ,φ) in the sample to the surface in direction to some point P:(x0,y0,z0).
    Note: M is given in spherical coordinates and P in Cartesian coordinates.

    The geometry of the sample is so that Oy is the symmetry axis (along the droplet stream direction),
    Oz is along the photon beam and Ox is the remaining axis in the orthogonal basis.
    θ is taken between Oz and the project of M onto the plane xOz, φ between M and Oy, and r = √(x^2+y^2+z^2)
    μ0 is the radius of the cylinder
"""
function d_sphere_P(r::Array{Cdouble,1},φ::Array{Cdouble,1},θ::Array{Cdouble,1},x0::Cdouble,y0::Cdouble,z0::Cdouble,μ0::Cdouble)
    dp = zeros(Cdouble,length(r),length(φ),length(θ)); # create matrix of disctances
    if (sqrt(x0^2+y0^2)>maximum(r))
        for i in 1:length(r)
            # only if the emission is coming from inside the sample (for sharp edge)
            if (r[i]<μ0) 
                for j in 1:length(φ)
                    for k in 1:length(θ)
                        # check that the coordinates are in the adequate ranges
                        if (φ[j]<0.0)
                            φj = -φ[j]
                            θk = θ[k] + π
                        elseif (φ[j]>π)
                            φj = 2π - φ[j]
                            θk = θ[k] + π
                        else
                            φj = φ[j]
                            θk = θ[k]
                        end
                        θk = θk %2π
                        ## M coordinates, spherical --> cartesian
                        xm = r[i]*cos(θk)*sin(φj);
                        ym = r[i]*sin(θk)*sin(φj);
                        zm = r[i]*cos(φj);

                        ## compute all cos/sin
                        cosω = (x0 - xm) / sqrt((x0-xm)^2 + (y0-ym)^2); # azimuthal angle of the direction MP
                        sinω = (y0 - ym) / sqrt((x0-xm)^2 + (y0-ym)^2);

                        # debugging
                        sinβ= (z0 - zm) / sqrt((x0-xm)^2 + (y0-ym)^2 + (z0-zm)^2); # polar angle of the direction MP is π/2-β
                        cosβ = sqrt((x0-xm)^2 + (y0-ym)^2) / sqrt((x0-xm)^2 + (y0-ym)^2 + (z0-zm)^2);

                        # cosα = sin(φ[j]) * cosβ * (cos(θ[k]) * cosω + sin(θ[k]) * sinω) + cos(φ[j]) * sinβ;
                        cosα = sin(φj) * cosβ * (cos(θk) * cosω + sin(θk) * sinω) + cos(φj) * sinβ;
                        
                        #if ((μ0^2)>=(r[i]^2 * (1 - cosα^2)))
                        dp[i,j,k] = -r[i] * cosα + sqrt(μ0^2 - r[i]^2 * (1 - cosα^2)); # distance from intersection to P
                        #end
                    end
                end
            else
                # no electron emitted from the shadow of the sample will reach the kinetic energy analyzer: set the distance to μ0 
                for j in 1:length(φ)
                    for k in 1:length(θ)
                        # check that the coordinates are in the adequate ranges
                        if (φ[j]<0.0)
                            φj = -φ[j]
                            θk = θ[k] + π
                        elseif (φ[j]>π)
                            φj = 2π - φ[j]
                            θk = θ[k] + π
                        else
                            φj = φ[j]
                            θk = θ[k]
                        end
                        θk = θk %2π
                        ## M coordinates, spherical --> cartesian
                        xm = r[i]*cos(θk)*sin(φj);
                        ym = r[i]*sin(θk)*sin(φj);
                        zm = r[i]*cos(φj);
                        
                        cosθ = (xm*x0 + ym*y0 + zm*z0)/(sqrt(xm^2+ym^2+zm^2)*sqrt(x0^2+y0^2+z0^2))
                        if cosθ<0.0
                            dp[i,j,k] = μ0
                        end
                    end
                end
            end
        end
    else
        throw("Not a good place for the analyzer: in the middle of the stream of droplets")
    end
    dp[dp.<0.0] .= 0.0;
    dp
end
function d_sphere_P(r::Array{Cdouble,1},φ::Array{Cdouble,1},θ::Cdouble,x0::Cdouble,y0::Cdouble,z0::Cdouble,μ0::Cdouble)
    dp = zeros(Cdouble,length(r),length(φ)); # create matrix of disctances
    if (sqrt(x0^2+y0^2)>maximum(r))
        ## M coordinates, spherical --> cartesian
        xm = cos(θ)*r*sin.(φ')
        ym = sin(θ)*r*sin.(φ')
        zm = r*cos.(φ')
        ## compute all cos/sin
        cosω = (x0 .- xm) ./ sqrt.((x0.-xm).^2 .+ (y0.-ym).^2); # azimuthal angle of the direction MP
        sinω = (y0 .- ym) ./ sqrt.((x0.-xm).^2 .+ (y0.-ym).^2)

        # debugging
        sinβ = (z0 .- zm) ./ sqrt.((x0.-xm).^2 .+ (y0.-ym).^2 .+ (z0.-zm).^2); # polar angle of the direction MP is π/2-β
        cosβ = sqrt.((x0.-xm).^2 .+ (y0.-ym).^2) ./ sqrt.((x0.-xm).^2 .+ (y0.-ym).^2 .+ (z0.-zm).^2);
        
        cosα =  cosβ .* (cos(θ)*cosω+sin(θ)*sinω) .* sin.(φ')  +  sinβ .* cos.(φ') ;

        A = -r.*cosα
        B = μ0^2 .- r.^2 .* (1.0 .- cosα.^2)
        A[r.>μ0,:] .= 0.0
        B[r.>μ0,:] .= 0.0
        dp = A+sqrt.(B);
        cosθ = (xm*x0 .+ ym*y0 .+ zm*z0)./(sqrt.(xm.^2+ym.^2+zm.^2)*sqrt(x0^2+y0^2+z0^2))
        dp[cosθ.<0.0] .= μ0
    else
        throw("Not a good place for the analyzer: in the middle of the stream of droplets")
    end
    dp[dp.<0.0] .= 0.0;
    dp
end



##
## discretized models
##
"""
    plane_gain_H(x::Array{Cdouble,1},y::Array{Cdouble,1},z::Array{Cdouble,1},x0::Cdouble,y0::Cdouble,z0::Cdouble,λe::Cdouble)

    Compute the volume integrales (exact for piecewise linear gain)

    H_{n,j,k} = ∭ e_n(z)e_j(x)e_k(y) e^{-d_P(M)/λe} dzdxdy

    The arrays x, y and z are the discretization subdivisions
    P:(x0,y0,z0) is the point in Cartesian coordinates used for the computation of the distance d_P
    λe is the attenuation length in the Beer-Lambert model

    Note: does not include the profile of the beam light, but it can be add easily if known
"""
function plane_gain_H(x::Array{Cdouble,1},y::Array{Cdouble,1},z::Array{Cdouble,1},x0::Cdouble,y0::Cdouble,z0::Cdouble,λe::Cdouble)
    # compute the elementary integrals
    Azn = 0.5*[z[2]-z[1]; z[3:end]-z[1:end-2]; z[end]-z[end-1]]; # dz
    Axj = 0.5*[x[2]-x[1]; x[3:end]-x[1:end-2]; x[end]-x[end-1]]; # dx
    Ayk = 0.5*[y[2]-y[1]; y[3:end]-y[1:end-2]; y[end]-y[end-1]]; # dy

    #compute the model
    H_zxy = zeros(Cdouble,length(z),length(x),length(y));

    for j in 1:length(x)
        for k in 1:length(y)
            H_zxy[:,j,k] = exp.(-d_plane_P(x[j],y[k],z,x0,y0,z0)/λe)
        end
    end

    H_z  = zeros(Cdouble,length(z));
    for n in 1:length(z)
        H_z[n]  = Azn[n]*Axj'*H_zxy[n,:,:]*Ayk
    end

    H_z,H_zxy,Azn,Axj,Ayk
end

"""
    finger_gain_H(x::Cdouble,y::Cdouble,z::Array{Cdouble,1},x0::Cdouble,y0::Cdouble,z0::Cdouble,λe::Cdouble)

    Compute the volume integrales (exact for piecewise linear gain)

    H_{n} = ∫ e_n(z) e^{z/λe} dz

    The array z is the discretization subdivisions, and (x,y) is the coordinate at the sample's surface
    where the model is evaluated.
    P:(x0,y0,z0) is the point in Cartesian coordinates used for the computation of the distance d_P
    λe is the attenuation length in the Beer-Lambert model

    Note: this model should be used only if not enough information is known about the experiment's geometry. It assumes
    that the whole signal is coming from one line orthogonal to the sample's surface, thus ignoring the variation in
    the distance the electron has to travel though the sample due to the surfacic spread of the illumination,
    making this model less surface sensitive.
"""
function finger_gain_H(x::Cdouble,y::Cdouble,z::Array{Cdouble,1},x0::Cdouble,y0::Cdouble,z0::Cdouble,λe::Cdouble)
    # compute the elementary integrals
    Azn = 0.5*[z[2]-z[1]; z[3:end]-z[1:end-2]; z[end]-z[end-1]]; # dz

    #compute the model
    H_z  = Azn.*exp.(-d_plane_P(x,y,z,x0,y0,z0)/λe);
    H_z,Azn
end

"""
    cylinder_gain_H(r::Array{Cdouble,1},θ::Array{Cdouble,1},y::Array{Cdouble,1},x0::Cdouble,y0::Cdouble,z0::Cdouble,μ0::Cdouble,λe::Cdouble)

    Compute the volume integrales (exact for piecewise linear gain)

    H_{n,j,k} = ∭ e_n(r)e_j(θ)e_k(y) e^{-d_P(M)/λe} rdrdθdy

    The arrays r, θ and y are the discretization subdivisions
    P:(x0,y0,z0) is the point in Cartesian coordinates used for the computation of the distance d_P
    μ0 is the radius of the cylinder
    λe is the attenuation length in the Beer-Lambert model

    Note: does not include the profile of the beam light, 
    but it can be add easily if known, see alignmentParameter
"""
function cylinder_gain_H(r::Array{Cdouble,1},θ::Array{Cdouble,1},y::Array{Cdouble,1},x0::Cdouble,y0::Cdouble,z0::Cdouble,μ0::Cdouble,λe::Cdouble)
    # compute the elementary integrals
    Arn = 0.5*r.*[r[2]-r[1]; r[3:end]-r[1:end-2]; r[end]-r[end-1]]; # rdr
    Aθj = 0.5*[θ[2]-θ[1]; θ[3:end]-θ[1:end-2]; θ[end]-θ[end-1]];    # dθ
    Ayk = 0.5*[y[2]-y[1]; y[3:end]-y[1:end-2]; y[end]-y[end-1]];    # dy

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

"""
sphere_gain_H(r::Array{Cdouble,1},θ::Array{Cdouble,1},φ::Array{Cdouble,1},x0::Cdouble,y0::Cdouble,z0::Cdouble,μ0::Cdouble,λe::Cdouble)

    Compute the volume integrales (exact for piecewise linear gain)

    H_{n,j,k} = ∭ e_n(r)e_j(θ)e_k(φ) e^{-d_P(M)/λe} r^2sin(φ)drdφdθ

    The geometry of the sample is so that Oy is the symmetry axis (along the droplet stream direction),
    Oz is along the photon beam and Ox is the remaining axis in the orthogonal basis.

    The arrays r, φ and θ are the discretization subdivisions (r: radial distance, φ polar angle, and θ azimuthal angle).
    The polar axis is Oy, φ is the polar angle between M and Oy, θ is taken between Oz and the project of M onto the plane xOz, and r = √(x^2+y^2+z^2)

    P:(x0,y0,z0) is the point in Cartesian coordinates used for the computation of the distance d_P
    μ0 is the radius of the sphere
    λe is the attenuation length in the Beer-Lambert model

    Note: does not include the profile of the beam light, 
    but it can be add easily if known, see alignmentParameterSphere
"""
function sphere_gain_H(r::Array{Cdouble,1},φ::Array{Cdouble,1},θ::Array{Cdouble,1},x0::Cdouble,y0::Cdouble,z0::Cdouble,μ0::Cdouble,λe::Cdouble)
    # compute the elementary integrals
    Ar1 = (1.0/(r[2]-r[1]))*( (r[2]/3.0)*(r[2]^3-r[1]^3) - (1.0/4.0)*(r[2]^4-r[1]^4) );
    ArN = (1.0/(r[end]-r[end-1]))*( (1.0/4.0)*(r[end]^4-r[end-1]^4) - (r[end-1]/3.0)*(r[end]^3-r[end-1]^3));
    Arn = [Ar1; (1.0/12.0)*((r[3:end].^3-r[1:end-2].^3) .+ r[2:end-1].*(r[3:end].^2-r[1:end-2].^2) .+ (r[3:end]-r[1:end-2])); ArN]; # r^2dr

    Aφ1 = cos(φ[1]) - (1.0/(φ[2]-φ[1]))*(sin(φ[2]) - sin(φ[1]));
    AφJ = (1.0/(φ[end]-φ[end-1]))*(sin(φ[end]) - sin(φ[end-1])) - cos(φ[end]);
    Aφj = [Aφ1; ((sin.(φ[2:end-1])-sin.(φ[1:end-2]))./(φ[2:end-1]-φ[1:end-2]) - (sin.(φ[3:end])-sin.(φ[2:end-1]))./(φ[3:end]-φ[2:end-1])); AφJ];    # sinφdφ

    Aθk = 0.5*[θ[2]-θ[1]; θ[3:end]-θ[1:end-2]; θ[end]-θ[end-1]];    # dθ

    #compute the model
    H_rφθ = zeros(Cdouble,length(r),length(φ),length(θ));
    for k in 1:length(θ)
        H_rφθ[:,:,k] = exp.(-d_sphere_P(r,φ,θ[k],x0,y0,z0,μ0)/λe)
    end

    H_r  = zeros(Cdouble,length(r));
    for n in 1:length(r)
        H_r[n] = Arn[n]*Aφj'*H_rφθ[n,:,:]*Aθk
    end
    H_r,H_rφθ,Arn,Aφj,Aθk
end


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


# """
#     Ψ_lin_peak(wsGeom::cylinderGeom,wsAcq::XPSacq;κ_cs::Cdouble=0.0,κ_eal::Cdouble=0.0)

#     returns the measurement operator for a cylindrical geometry
#     - wsGeom: cylinderGeom
#     - wsAcq:  XPSacq
#     - κ_cs:   relative bias in the cross section values
#     - κ_eal:  relative bias in the attenuation length values
# """
# function Ψ_lin_peak(wsGeom::cylinderGeom,wsAcq::XPSacq;κ_cs::Cdouble=0.0,κ_eal::Cdouble=0.0)
#     Hr,_,_,_,_ = cylinder_gain_H(wsGeom.r,wsGeom.θ,wsGeom.y,wsGeom.x0,wsGeom.y0,wsGeom.z0,wsGeom.μ0,(1.0+κ_eal)*wsAcq.λe);
#     wsAcq.T*wsAcq.α*wsAcq.Fν*((1.0+κ_cs)*wsAcq.σν)*Hr'
# end

# """
#     Ψ_lin_peak_area(wsGeom::cylinderGeom,wsAcq::XPSacq;κ_cs::Cdouble=0.0,κ_eal::Cdouble=0.0)

#     returns the peak area operator for a cylindrical geometry
#     - wsGeom: cylinderGeom
#     - wsAcq:  XPSacq
#     - κ_cs:   relative bias in the cross section values
#     - κ_eal:  relative bias in the attenuation length values
# """
# function Ψ_lin_peak_area(wsGeom::cylinderGeom,wsAcq::XPSacq;κ_cs::Cdouble=0.0,κ_eal::Cdouble=0.0)
#     H_peak = Ψ_lin_peak(wsGeom,wsAcq;κ_cs=κ_cs,κ_eal=κ_eal)
#     AKe = 0.5*[wsAcq.Ke[2]-wsAcq.Ke[1]; wsAcq.Ke[3:end]-wsAcq.Ke[1:end-2]; wsAcq.Ke[end]-wsAcq.Ke[end-1]];    # dθ
#     AKe'*H_peak # A = dropdims(sum(H,dims=1),dims=1)*abs(peak.Ke[2]-peak.Ke[1])
# end


"""
    cov_H_cylinder()

    computes the covariance matrix of the geometrical structure...
    A = ∭ ρ(r,θ,y) e^{\\frac{d_P(r,θ,y)}{λ}} r dr dθ dy ≃ Hρ
    where H = [H_1 H_2 … H_Nr]† and 
    H_n(λ) = r_n ∑_j ∑_k e^{\\frac{d_P(r_n,θ_j,y_k)}{λ}} ∭ e_n(r) e_j(θ) e_k(y) dr dθ dy
    ΓH = cov(H) = \\mathbb{E} [H×H†] - \\mathbb{E}[H]×\\mathbb{E} [H†]
"""
function cov_H_cylinder(r::Array{Cdouble,1},θ::Array{Cdouble,1},y::Array{Cdouble,1},x0::Cdouble,y0::Cdouble,z0::Cdouble,μ0::Cdouble,λ::Array{Cdouble,1},Pλ::Array{Cdouble,1})
    # the distance for each point of the space discretization
    Nr = length(r);
    Nθ = length(θ);
    Ny = length(y);
    D = zeros(Cdouble,Nr,Nθ,Ny);
    for k in 1:Ny
        D[:,:,k] = d_cylinder_P(r,θ,y[k],x0,y0,z0,μ0);
    end
    Arn = 0.5*r.*[r[2]-r[1]; r[3:end]-r[1:end-2]; r[end]-r[end-1]]; # rdr
    Aθj = 0.5*[θ[2]-θ[1]; θ[3:end]-θ[1:end-2]; θ[end]-θ[end-1]];    # dθ
    Ayk = 0.5*[y[2]-y[1]; y[3:end]-y[1:end-2]; y[end]-y[end-1]];    # dy

    # attenuation length distribution
    Nλ = length(λ);
    Aλ = 0.5*[λ[2]-λ[1]; λ[3:end]-λ[1:end-2]; λ[end]-λ[end-1]];

    # compute the integration operator for the discretized attenuation space
    H = zeros(Cdouble,Nr,Nλ);
    Djk = Aθj*Ayk'; # integration over θ and y 
    for m in 1:Nr
        for s in 1:Nλ
            H[m,s] = Arn[m]*sum(Djk.*exp.(-D[m,:,:]/λ[s]))
        end
    end

    # mean operator 
    μH = H*(Pλ.*Aλ);

    # square
    HHt = zeros(Cdouble,Nr,Nr);
    for l in 1:Nr
        HHt[l,l] = (H[l,:].^2 .*Pλ)'*Aλ
        for m in l+1:Nr
            HHt[l,m] = (H[l,:].*H[m,:].*Pλ)'*Aλ
            HHt[m,l] = HHt[l,m]
        end
    end

    # return the covariance
    HHt - μH*μH', μH # , H
end

# TODO: compute the model's uncertainty (the ones due to errors in the parameters' values)
# # uniform distribution
# function Ψ_lin_peaks_mean_and_std(wsGeom::cylinderGeom,wsAcq::XPSacq;κ_cs::Cdouble=0.0,κ_eal::Cdouble=0.0)
#     N = length(Zi);
#     H_mean = Array{Cdouble}(undef,wsXPS.Nν*wsXPS.Nbe,N);
#     H_var  = Array{Cdouble}(undef,wsXPS.Nν*wsXPS.Nbe,N);
#     # go for the discrete operator
#     for j in 1:wsXPS.Nν
#        for b in 1:wsXPS.Nbe
#           λi_min = (1.0-κ_eal)*wsXPS.λe[(j-1)*wsXPS.Nbe+b];
#           λi_max = (1.0+κ_eal)*wsXPS.λe[(j-1)*wsXPS.Nbe+b];
#           H_mean[(j-1)*wsXPS.Nbe+b,1] = wsXPS.αT[j]*wsXPS.Fν[j]*wsXPS.σν[j,b]*f_ij_0(Zi[1],Zi[2],λi_min,λi_max,σ_z,z0,Nz,Nλ)
#           H_var[(j-1)*wsXPS.Nbe+b,1]  = (wsXPS.αT[j]*wsXPS.Fν[j]*wsXPS.σν[j,b])^2*(1.0+(κ_cs^2/3.0))*f_ij_0_sq(Zi[1],Zi[2],λi_min,λi_max,σ_z,z0,Nz,Nλ)
#           for i in 2:N-1
#              H_mean[(j-1)*wsXPS.Nbe+b,i] = wsXPS.αT[j]*wsXPS.Fν[j]*wsXPS.σν[j,b]*f_ij(Zi[i-1],Zi[i],Zi[i+1],λi_min,λi_max,σ_z,z0,Nz,Nλ)
#              H_var[(j-1)*wsXPS.Nbe+b,i]  = (wsXPS.αT[j]*wsXPS.Fν[j]*wsXPS.σν[j,b])^2*(1.0+(κ_cs^2/3.0))*f_ij_sq(Zi[i-1],Zi[i],Zi[i+1],λi_min,λi_max,σ_z,z0,Nz,Nλ)
#           end
#           H_mean[(j-1)*wsXPS.Nbe+b,end] = wsXPS.αT[j]*wsXPS.Fν[j]*wsXPS.σν[j,b]*f_ij_M(Zi[end-1],Zi[end],λi_min,λi_max,σ_z,z0,Nz,Nλ)
#           H_var[(j-1)*wsXPS.Nbe+b,end]  = (wsXPS.αT[j]*wsXPS.Fν[j]*wsXPS.σν[j,b])^2*(1.0+(κ_cs^2/3.0))*f_ij_M_sq(Zi[end-1],Zi[end],λi_min,λi_max,σ_z,z0,Nz,Nλ)
#        end
#     end
#     H_mean,sqrt.(H_var-H_mean.^2)
#  end


#  """
#     Ψ_lin_peaks(wsGeom::cylinderGeom,XPS_peak::Dict{Int64,XPSacq};κ_cs::Cdouble=0.0,κ_eal::Cdouble=0.0)

#     return a dictionary of measurement operators corresponding to the peaks listed in the argument

#     - wsGeom::cylinderGeom    
#     - XPS_peak: dictionary of peaks (different acquisition setups)
#     - κ_cs:   relative bias in the cross section values
#     - κ_eal:  relative bias in the attenuation length values
#  """
# function Ψ_lin_peaks(wsGeom::cylinderGeom,XPS_peak::Dict{Int64,XPSacq};κ_cs::Cdouble=0.0,κ_eal::Cdouble=0.0)
#     H_dict     = Dict{Int64,Array{Cdouble,2}}();
    
#     for (i,peak) in XPS_peak
#        # create the measurement model
#        setindex!(H_dict,Ψ_lin_peak(wsGeom,peak;κ_cs=κ_cs,κ_eal=κ_eal),i);
#     end
#     H_dict
# end

