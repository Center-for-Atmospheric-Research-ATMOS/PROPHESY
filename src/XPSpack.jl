# push!(LOAD_PATH,"/path/to/package/XPSinv.jl/src/")
# sudo nano /opt/julias/julia-1.3.1/etc/julia/startup.jl


module XPSpack

using Statistics
using LinearAlgebra # for the definition of the I matrix in XPSutils... could be avoided, but the package is very common
using Interpolations # used in XPSmeas.jl
using Printf
using utilsFun

# export the functions for the end user
export λe, σ_cs_orb, σ_exp, Ψ, Ψ_lin, Ψ_lin_0, Ψ_lin_M # electron attenuation length, corss section, measurement operators
export Ψ_peaks, Ψ_lin_peaks                  # should become the main functions for the end user

# uncertainty: model uncertainty due to parameter uncertainty
export Ψij_distribution!, mean_and_std_map   # variability in the measurement operator due to variability in λe and σ_exp
export mean_and_std_map_peaks, mean_and_std_map_lin_peaks
export Ψ_lin_peaks_mean_and_std
export Ψij_distribution_lin!, Ψij_distribution_lin_0!, Ψij_distribution_lin_M!, mean_and_std_map_lin # more variability estimation

# other uncertainty
export δintegrals                            # for the estimation of the error due to the discretization of the concentraion profiles
export ΔΨλ, ΔΨσ                              # uncertainty: first order model

# other functions (parametric total concentration, riemann numerical integration and basis functions)
export ρ_tot_int, riemann, e_k, e_0, e_M     # other

# data enhancement: baseline correction and peak fitting (requires sampling of the spectra)
export baseline_removal, EM_peaks, cross_section_spread_function # in XPSutils.jl

# objects modelling experiment and device
export XPSexp, XPSdevice, XPSsetup

# include the implementation of the exported functions and objects
include("XPSexp.jl")   # implements XPSexp object
include("XPSmeas.jl")  # implement most function exported so far

using NewtonMethod
include("XPSutils.jl") # common algorithms used for data processing #TODO: move to XPSinv package

end # module
