# smooth_edge.jl

"""
    distance till the edge of the smooth sample

    τ_root(α::Cdouble,γ::Cdouble,x0::Cdouble,z0::Cdouble,zm::Cdouble,xm::Cdouble,R2::Cdouble)

    input: 

      - α and γ: angular direction along the path between M=(xm,ym,zm) and P=(x0,y0,z0)
      - x0 and z0: cartesian coordinate of the analyzer in the plan zOx
      - xm and zm: cartesian coordinate of M in the plan zOx
      - R2: difference of square radial distance between the edge of the sample and M, R2 = (μ0+κΔr)^2 - rm^2 ⩾0

    output: 

      - return the distance between the current location (zm,xm) and the end of the edge of the sample in the direction of the aperture of the analyzer in (z0,x0)
"""
function τ_root(α::Cdouble,γ::Cdouble,x0::Cdouble,z0::Cdouble,zm::Cdouble,xm::Cdouble,R2::Cdouble)
    0.0
end

"""
    smooth edge model for the total concentration

    ρ_T_smooth(r::Array{Cdouble,1},μ0::Cdouble,Δr::Cdouble;ρ0::Cdouble=1.0)

    input:

      - r: radial distance to the center 
      - μ0: radius of the sample
      - Δr: transition length
      - ρ0: total bulk concentration

    output: 

      - total concentration 
"""
function ρ_T_smooth(r::Array{Cdouble,1},μ0::Cdouble,Δr::Cdouble;ρ0::Cdouble=1.0)
    0.0
end

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
    0.0
end