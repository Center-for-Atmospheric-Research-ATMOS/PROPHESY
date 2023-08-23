# push!(LOAD_PATH,"/path/to/package/XPSpack/")
# sudo nano /opt/julias/julia-1.3.1/etc/julia/startup.jl
# or nano ~/.julia/julia-1.8.0/etc/julia/startup.jl
# or nano ~/.julia/config/startup.jl

"""
This is the [`XPSpack`](@ref), which contains
* [`XPSpack.σ_bg`](@ref)
* [`XPSpack.σ_bg_density`](@ref)
* [`XPSpack.fingerGeom`](@ref)
* [`XPSpack.planeGeom`](@ref)
* [`XPSpack.cylinderGeom`](@ref)
* [`XPSpack.d_plane_P`](@ref)
* [`XPSpack.d_cylinder_P`](@ref)
* [`XPSpack.d_cylinder_P_simple`](@ref)
* [`XPSpack.d_sphere_P`](@ref)
* [`XPSpack.plane_gain_H`](@ref)
* [`XPSpack.finger_gain_H`](@ref)
* [`XPSpack.cylinder_gain_H`](@ref)
* [`XPSpack.sphere_gain_H`](@ref)
* [`XPSpack.beamProfile`](@ref)
* [`XPSpack.alignmentParameter`](@ref)
* [`XPSpack.alignmentParameterSphere`](@ref)
* [`XPSpack.alignmentParameterSmoothEdge`](@ref)
* [`XPSpack.alignmentParameterSphereSmoothEdge`](@ref)
* [`XPSpack.cov_H_cylinder`](@ref)
* [`XPSpack.baseline_removal`](@ref)
* [`XPSpack.EM_peaks`](@ref)
* [`XPSpack.cross_section_spread_function`](@ref)
* [`XPSpack.cross_section_spread_function_sample`](@ref)
* [`XPSpack.noiseAndParameterEstimation`](@ref)
* [`XPSpack.cylinder_gain_smooth_H`](@ref)
* [`XPSpack.cov_H_cylinder_smooth`](@ref)
* [`XPSpack.sphere_gain_smooth_H`](@ref)
* [`XPSpack.cov_H_sphere_smooth`](@ref)
"""
module XPSpack

using Statistics
using LinearAlgebra # for the definition of the I matrix in XPSutils
using Interpolations # used in penetration_depth.jl 
using Printf
using NMOpt    # package at (https://github.com/matthewozon/NMOpt)    # used in XPSutils.jl
using MINOTAUR # package at (https://github.com/matthewozon/MINOTAUR) # used in smooth_edge.jl for the computation of integrals (quadrature of composed function)

# penetration depth and cross-section values: not sure it's really useful in this module since the values are supposed to be known to create the models
export λe, λe_exp, σ_cs_orb, σ_bg_density, σ_bg_lin_density, σ_bg
export σ_C1s_exp,σ_O1s_exp,σ_S2p_exp

# data enhancement: baseline correction and peak fitting (requires sampling of the spectra)
export baseline_removal, EM_peaks, cross_section_spread_function, cross_section_spread_function_sample # in XPSutils.jl
export logistic

# possible geometry of the sample
export fingerGeom, planeGeom, cylinderGeom
# distance and geometry factors
export d_plane_P, d_sphere_P, plane_gain_H, finger_gain_H
export cylinder_gain_H, d_cylinder_P, d_cylinder_P_simple
export sphere_gain_H,d_sphere_P
export beamProfile

# smooth edge model
export cylinder_gain_smooth_H, sphere_gain_smooth_H

# uncertainty
export cov_H_cylinder, cov_H_cylinder_smooth, cov_H_sphere_smooth

# function for noise estimation
export noiseAndParameterEstimation                                      # estimation based on noise model
export alignmentParameter, alignmentParameterSphere                     # definition using the sharp edge model
export alignmentParameterSmoothEdge, alignmentParameterSphereSmoothEdge # definition using the smooth edge model

# include the implementation of the exported functions and objects
include("penetration_depth.jl")
include("cross_section.jl")
include("XPSmeas.jl")  # implement most function exported so far
include("smooth_edge.jl")
include("XPSalignment.jl") # alignment parameter calculations
include("XPSutils.jl") # common algorithms used for data processing
include("geom_struct.jl") # implements fingerGeom, planeGeom, cylinderGeom (should be deprecated)

end # module
