#------------------------------------------------------------------------------
#
# This file is part of the XPSpack module which is licensed under CC-BY 4.0 license.
#
# Copyright (C) 2022,  Matthew Ozon, Konstantin Tumashevich and Nønne L. Prisle.
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

    The geometry of the sample is so that Oz is the symmetry axis (along the droplet stream direction),
    Ox is along the photon beam and Oy is the remaining axis in the orthogonal basis.
    θ is taken between Ox and the projection of M onto the plane xOy, φ between M and Oz, and r = √(x^2+y^2+z^2)
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

    The geometry of the sample is so that Oz is the symmetry axis (along the droplet stream direction),
    Ox is along the photon beam and Oy is the remaining axis in the orthogonal basis.

    The arrays r, φ and θ are the discretization subdivisions (r: radial distance, φ polar angle, and θ azimuthal angle).
    The polar axis is Oz, φ is the polar angle between M and Oz, θ is taken between Ox and the projection of M onto the plane xOy, and r = √(x^2+y^2+z^2)

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
