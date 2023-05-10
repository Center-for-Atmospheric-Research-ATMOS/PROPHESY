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
* [`XPSpack.cov_H_cylinder`](@ref)
* [`XPSpack.baseline_removal`](@ref)
* [`XPSpack.EM_peaks`](@ref)
* [`XPSpack.cross_section_spread_function`](@ref)
* [`XPSpack.cross_section_spread_function_sample`](@ref)
* [`XPSpack.noiseAndParameterEstimation`](@ref)
* [`XPSpack.τ_root`](@ref)
* [`XPSpack.ρ_T_smooth`](@ref)
* [`XPSpack.cylinder_gain_smooth_H`](@ref)
"""
module XPSpack

using Statistics
using LinearAlgebra # for the definition of the I matrix in XPSutils
using Interpolations # used in penetration_depth.jl 
using Printf
using NMOpt # unregistered package at (https://github.com/matthewozon/NMOpt)  # used in XPSutils.jl

# penetration depth and cross-section values: not sure it's really useful in this module since the values are supposed to be known to create the models
export λe, λe_exp, σ_cs_orb, σ_bg_density, σ_bg_lin_density, σ_bg
export σ_C1s_exp,σ_O1s_exp,σ_S2p_exp

# data enhancement: baseline correction and peak fitting (requires sampling of the spectra)
export baseline_removal, EM_peaks, cross_section_spread_function, cross_section_spread_function_sample # in XPSutils.jl
export logistic

# possible geometry of the sample
export fingerGeom, planeGeom, cylinderGeom # not used: should be deprecated
# distance and geometry factors
export d_plane_P, d_sphere_P, plane_gain_H, finger_gain_H
export cylinder_gain_H, d_cylinder_P, d_cylinder_P_simple
export cov_H_cylinder
export sphere_gain_H,d_sphere_P
export beamProfile

# smooth edge model
export cylinder_gain_smooth_H, ρ_T_smooth, τ_root

# function for noise estimation
export alignmentParameter, alignmentParameterSphere, noiseAndParameterEstimation

# include the implementation of the exported functions and objects
include("penetration_depth.jl")
include("cross_section.jl")
include("XPSmeas.jl") # measurement models
include("smooth_edge.jl")
include("XPSalignment.jl") # alignment parameter calculations
include("XPSutils.jl") # common algorithms used for data processin
include("geom_struct.jl") # implements fingerGeom, planeGeom, cylinderGeom (should be deprecated)

end # module
