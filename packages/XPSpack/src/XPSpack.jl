# push!(LOAD_PATH,"/path/to/package/XPSpack/")
# sudo nano /opt/julias/julia-1.3.1/etc/julia/startup.jl
# or nano ~/.julia/julia-1.8.0/etc/julia/startup.jl
# or nano ~/.julia/config/startup.jl


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
export fingerGeom, planeGeom, cylinderGeom
# distance and geometry factors
export d_plane_P, d_sphere_P, plane_gain_H, finger_gain_H
export cylinder_gain_H,alignmentParameter, d_cylinder_P, d_cylinder_P_simple
export cov_H_cylinder
export sphere_gain_H,d_sphere_P,alignmentParameterSphere
export beamProfile

# deprecation before merging beta to main
# objects modelling experiment and device
# export XPSacq # acquisition parameters
# export Ψ_lin_peaks # implemented in XPSmeas.jl # used in main branch in: data_gen_exp_5.jl, SDS_water.jl, data_generation_exp.jl 

# function for noise estimation
export noiseAndParameterEstimation

# include the implementation of the exported functions and objects
include("penetration_depth.jl")
include("cross_section.jl")
include("XPSmeas.jl")  # implement most function exported so far
include("XPSutils.jl") # common algorithms used for data processing


end # module
