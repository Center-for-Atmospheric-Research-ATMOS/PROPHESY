# push!(LOAD_PATH,"/path/to/package/XPSinv.jl/src/")
# sudo nano /opt/julias/julia-1.3.1/etc/julia/startup.jl


module XPSpack

using Statistics
using LinearAlgebra # for the definition of the I matrix in XPSutils... could be avoided, but the package is very common
using Interpolations # used in XPSmeas.jl
using Printf
using utilsFun
using NewtonMethod  # used in XPSutils.jl

# penetration depth and cross-section values: not sure it's really useful in this module since the values are supposed to be known to create the models
export λe, σ_cs_orb 
include("XPSpack/penetration_depth.jl")
include("XPSpack/cross_section.jl")


# data enhancement: baseline correction and peak fitting (requires sampling of the spectra)
export baseline_removal, EM_peaks, cross_section_spread_function, cross_section_spread_function_sample # in XPSutils.jl


# possible geometry of the sample
export fingerGeom, planeGeom, cylinderGeom
# distance and geometry factors
export d_plane_P, d_cylinder_P, d_cylinder_P_simple, d_sphere_P, plane_gain_H, finger_gain_H, cylinder_gain_H
export cov_H_cylinder

# objects modelling experiment and device
export XPSacq # acquisition parameters
export Ψ_lin_peaks # implemented in XPSmeas.jl

# include the implementation of the exported functions and objects
include("XPSpack/XPSmeas.jl")  # implement most function exported so far


include("XPSpack/XPSutils.jl") # common algorithms used for data processing

# MAYBE: move to another package/module
export samplePosterior, acceptSample, transmissionMechanism, smoothnessCovariance, corrCovariance
export samplePosteriorModelMargin, acceptSampleModelMargin # marginalization over the measurement operator space (or some approximation of it)
export samplePosteriorEntropy, acceptSampleEntropy # does ot serve much purpose (except for showing that it's not gonna be used)
export samplePosteriorMargin, acceptSampleMargin
include("XPSsampling.jl")

end # module
