"""
    smooth edge model attenuation gain

    cylinder_gain_smooth_H(r::Array{Cdouble,1},θ::Array{Cdouble,1},y::Array{Cdouble,1},x0::Cdouble,y0::Cdouble,z0::Cdouble,μ0::Cdouble,Δr::Cdouble,λe::Cdouble;κ::Cdouble=3.0,Nτ::Int64=20)

    input:

      - r,θ,y: cylindrical coordinates of the discretization of the cylinder (symmetry axis Oy)
      - x0,y0,z0: cartesian coordinate of the analyzer aperture
      - μ0: radius of the sample (equivalent radius of the sharp edge model)
      - Δr: transition distance (at r = μ0+Δr, the total concentration is ρ0/(1+e)≃0.269ρ0, and at r = μ0-Δr, the concentration is ρ0×e/(1+e)≃0.731ρ0)
      - λe: attenuation length
    
    optional input:

      - κ:  the maximum radial distance of the attenuation integral is set at r=μ0+κΔr
      - Nτ: number of discretization points for the attenuation integral

    output:

      - H_r: integrated gain (depends only the radial distance)
      - H_rθy: 3D array containing the attenuation gain evaluated for the discretization nodes (r,θ,y)
      - Arn,Aθj,Ayk: discretization of the volume integrals
""" 
function cylinder_gain_smooth_H(r::Array{Cdouble,1},θ::Array{Cdouble,1},y::Array{Cdouble,1},x0::Cdouble,y0::Cdouble,z0::Cdouble,μ0::Cdouble,Δr::Cdouble,λe::Cdouble;κ::Cdouble=5.0,Nτ::Int64=20)
    # geometry factor gain for each discretization point
    H_rθy = exp.(-quadrature_fg_cylinder(r,θ,y,x0,y0,z0,μ0,Δr;κ=κ,Nτ=Nτ)/λe);

    # compute the elementary integrals
    Arn = 0.5*r.*[r[2]-r[1]; r[3:end]-r[1:end-2]; r[end]-r[end-1]]; # rdr
    Aθj = 0.5*[θ[2]-θ[1]; θ[3:end]-θ[1:end-2]; θ[end]-θ[end-1]];    # dθ
    Ayk = 0.5*[y[2]-y[1]; y[3:end]-y[1:end-2]; y[end]-y[end-1]];    # dy

    H_r  = zeros(Cdouble,length(r));
    for n in 1:length(r)
        H_r[n] = Arn[n]*Aθj'*H_rθy[n,:,:]*Ayk
    end

    # return 
    H_r,H_rθy,Arn,Aθj,Ayk
end


"""
    cov_H_cylinder_smooth(r::Array{Cdouble,1},θ::Array{Cdouble,1},y::Array{Cdouble,1},x0::Cdouble,y0::Cdouble,z0::Cdouble,μ0::Cdouble,λ::Array{Cdouble,1},Pλ::Array{Cdouble,1},Δr::Cdouble)

    computes the covariance matrix of the geometrical structure...
    A = ∭ ρ(r,θ,y) e^{ ∫_{0}^{̄τ} \\frac{ρ(M_s(τ))}{ρ0*λ} dτ } r dr dθ dy ≃ Hρ
    where H = [H_1 H_2 … H_Nr]† and 
    H_n(λ) = r_n ∑_j ∑_k e^{ ∫_{0}^{̄τ} \\frac{ρ(M_s(τ;r_n,θ_j,y_k))}{ρ0*λ} dτ } ∭ e_n(r) e_j(θ) e_k(y) dr dθ dy
    ΓH = cov(H) = \\mathbb{E} [H×H†] - \\mathbb{E}[H]×\\mathbb{E} [H†]
"""
function cov_H_cylinder_smooth(r::Array{Cdouble,1},θ::Array{Cdouble,1},y::Array{Cdouble,1},x0::Cdouble,y0::Cdouble,z0::Cdouble,μ0::Cdouble,Δr::Cdouble,λ::Array{Cdouble,1},Pλ::Array{Cdouble,1};κ::Cdouble=5.0,Nτ::Int64=20)
    # attenuation integral (or weighted distance) for each point of the space discretization
    D = quadrature_fg_cylinder(r,θ,y,x0,y0,z0,μ0,Δr;κ=κ,Nτ=Nτ)

    # compute the elementary integrals
    Arn = 0.5*r.*[r[2]-r[1]; r[3:end]-r[1:end-2]; r[end]-r[end-1]]; # rdr
    Aθj = 0.5*[θ[2]-θ[1]; θ[3:end]-θ[1:end-2]; θ[end]-θ[end-1]];    # dθ
    Ayk = 0.5*[y[2]-y[1]; y[3:end]-y[1:end-2]; y[end]-y[end-1]];    # dy

    # attenuation length distribution
    Nλ = length(λ);
    Aλ = 0.5*[λ[2]-λ[1]; λ[3:end]-λ[1:end-2]; λ[end]-λ[end-1]];

    # compute the integration operator for the discretized attenuation space
    Nr = length(r);
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
    HHt - μH*μH', μH 
end


"""
    smooth edge model attenuation gain

    sphere_gain_smooth_H(r::Array{Cdouble,1},φ::Array{Cdouble,1},θ::Array{Cdouble,1},x0::Cdouble,y0::Cdouble,z0::Cdouble,μ0::Cdouble,Δr::Cdouble,λe::Cdouble;κ::Cdouble=5.0,Nτ::Int64=20)

    input:

      - r,φ,θ: spherical coordinates of the discretization of the sphere (polar axis Oz)
      - x0,y0,z0: cartesian coordinate of the analyzer aperture
      - μ0: radius of the sample (equivalent radius of the sharp edge model)
      - Δr: transition distance (at r = μ0+Δr, the total concentration is ρ0/(1+e)≃0.269ρ0, and at r = μ0-Δr, the concentration is ρ0×e/(1+e)≃0.731ρ0)
      - λe: attenuation length
    
    optional input:

      - κ:  the maximum radial distance of the attenuation integral is set at r=μ0+κΔr
      - Nτ: number of discretization points for the attenuation integral

    output:

      - H_r: integrated gain (depends only the radial distance)
      - H_rφθ: 3D array containing the attenuation gain evaluated for the discretization nodes (r,φ,θ)
      - Arn,Aφj,Aθk: discretization of the volume integrals
""" 
function sphere_gain_smooth_H(r::Array{Cdouble,1},φ::Array{Cdouble,1},θ::Array{Cdouble,1},x0::Cdouble,y0::Cdouble,z0::Cdouble,μ0::Cdouble,Δr::Cdouble,λe::Cdouble;κ::Cdouble=5.0,Nτ::Int64=20)
    # geometry factor gain for each discretization point
    H_rφθ = exp.(-quadrature_fg_sphere(r,φ,θ,x0,y0,z0,μ0,Δr;κ=κ,Nτ=Nτ)/λe);

    # compute the elementary integrals
    Ar1 = (1.0/(r[2]-r[1]))*( (r[2]/3.0)*(r[2]^3-r[1]^3) - (1.0/4.0)*(r[2]^4-r[1]^4) );
    ArN = (1.0/(r[end]-r[end-1]))*( (1.0/4.0)*(r[end]^4-r[end-1]^4) - (r[end-1]/3.0)*(r[end]^3-r[end-1]^3));
    Arn = [Ar1; (1.0/12.0)*((r[3:end].^3-r[1:end-2].^3) .+ r[2:end-1].*(r[3:end].^2-r[1:end-2].^2) .+ (r[3:end]-r[1:end-2])); ArN]; # r^2dr

    Aφ1 = cos(φ[1]) - (1.0/(φ[2]-φ[1]))*(sin(φ[2]) - sin(φ[1]));
    AφJ = (1.0/(φ[end]-φ[end-1]))*(sin(φ[end]) - sin(φ[end-1])) - cos(φ[end]);
    Aφj = [Aφ1; ((sin.(φ[2:end-1])-sin.(φ[1:end-2]))./(φ[2:end-1]-φ[1:end-2]) - (sin.(φ[3:end])-sin.(φ[2:end-1]))./(φ[3:end]-φ[2:end-1])); AφJ];    # sinφdφ

    Aθk = 0.5*[θ[2]-θ[1]; θ[3:end]-θ[1:end-2]; θ[end]-θ[end-1]];    # dθ


    H_r  = zeros(Cdouble,length(r));
    for n in 1:length(r)
        H_r[n] = Arn[n]*Aφj'*H_rφθ[n,:,:]*Aθk
    end

    # return 
    H_r,H_rφθ,Arn,Aφj,Aθk
end

"""
    cov_H_sphere_smooth(r::Array{Cdouble,1},φ::Array{Cdouble,1},θ::Array{Cdouble,1},x0::Cdouble,y0::Cdouble,z0::Cdouble,μ0::Cdouble,Δr::Cdouble,λ::Array{Cdouble,1},Pλ::Array{Cdouble,1};κ::Cdouble=5.0,Nτ::Int64=20)

    computes the covariance matrix of the geometrical structure...
    A = ∭ ρ(r,φ,θ) e^{ ∫_{0}^{̄τ} \\frac{ρ(M_s(τ))}{ρ0*λ} dτ } r^2 sin φ dr dφ dθ  ≃ Hρ
    where H = [H_1 H_2 … H_Nr]† and 
    H_n(λ) = r_n^2 ∑_j sin φ_j  ∑_k e^{ ∫_{0}^{̄τ} \\frac{ρ(M_s(τ;r_n,φ_j,θ_k))}{ρ0*λ} dτ } ∭ e_n(r) e_j(φ) e_k(θ) dr dθ dy
    ΓH = cov(H) = \\mathbb{E} [H×H†] - \\mathbb{E}[H]×\\mathbb{E} [H†]
"""
function cov_H_sphere_smooth(r::Array{Cdouble,1},φ::Array{Cdouble,1},θ::Array{Cdouble,1},x0::Cdouble,y0::Cdouble,z0::Cdouble,μ0::Cdouble,Δr::Cdouble,λ::Array{Cdouble,1},Pλ::Array{Cdouble,1};κ::Cdouble=5.0,Nτ::Int64=20)
    # attenuation integral (or weighted distance) for each point of the space discretization
    D = quadrature_fg_sphere(r,φ,θ,x0,y0,z0,μ0,Δr;κ=κ,Nτ=Nτ);

    # compute the elementary integrals
    Ar1 = (1.0/(r[2]-r[1]))*( (r[2]/3.0)*(r[2]^3-r[1]^3) - (1.0/4.0)*(r[2]^4-r[1]^4) );
    ArN = (1.0/(r[end]-r[end-1]))*( (1.0/4.0)*(r[end]^4-r[end-1]^4) - (r[end-1]/3.0)*(r[end]^3-r[end-1]^3));
    Arn = [Ar1; (1.0/12.0)*((r[3:end].^3-r[1:end-2].^3) .+ r[2:end-1].*(r[3:end].^2-r[1:end-2].^2) .+ (r[3:end]-r[1:end-2])); ArN]; # r^2dr

    Aφ1 = cos(φ[1]) - (1.0/(φ[2]-φ[1]))*(sin(φ[2]) - sin(φ[1]));
    AφJ = (1.0/(φ[end]-φ[end-1]))*(sin(φ[end]) - sin(φ[end-1])) - cos(φ[end]);
    Aφj = [Aφ1; ((sin.(φ[2:end-1])-sin.(φ[1:end-2]))./(φ[2:end-1]-φ[1:end-2]) - (sin.(φ[3:end])-sin.(φ[2:end-1]))./(φ[3:end]-φ[2:end-1])); AφJ];    # sinφdφ

    Aθk = 0.5*[θ[2]-θ[1]; θ[3:end]-θ[1:end-2]; θ[end]-θ[end-1]];    # dθ

    # attenuation length distribution
    Nλ = length(λ);
    Aλ = 0.5*[λ[2]-λ[1]; λ[3:end]-λ[1:end-2]; λ[end]-λ[end-1]];

    # compute the integration operator for the discretized attenuation space
    Nr = length(r);
    H = zeros(Cdouble,Nr,Nλ);
    Djk = Aφj*Aθk'; # integration over φ and θ
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
    HHt - μH*μH', μH 
end

