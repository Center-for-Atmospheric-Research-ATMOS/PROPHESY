#
# XPSexp.jl --
#
# XPSexp.jl is part of the XPSmeas.jl module.
# It is an object containing the parameters used for the modellization of and
# XPS measurement: photon energy, photon flux, device characteristics, angle
# kinetic energies, cross sections and electron attenuation lengths
#
#------------------------------------------------------------------------------
#
# This file is part of the XPSexp module which is licensed under the MIT "Expat" License:
#
# Copyright (C) 2021-2022,  Matthew Ozon.
#
#-----------------------------------------------------------------------------


# #NOTE: I don't include the concentration profiles in this structure because the depth discretization is not a parameter of the experiment, it's only a modelling artifact
# mutable struct XPSexp

#     # setup
#     ħν::Cdouble                    # photon energy
#     Fν::Cdouble                    # photon flux
#     αT::Cdouble                    # device and experiment parameter (product of α and T)
#     θ::Cdouble                     # measurement angle (only one per experiment) (θ=0, sensor aligned with the polarization angle of the photon beam) (for now, just a tag)

#     # electron analyzer range
#     Nke::Int64                     # number of sampled kinetic energies
#     Ke::Array{Cdouble,1}           # kinetic energies

#     # the physics of the experiment (cross section and attenuation length)
#     σ_β::Array{Cdouble,1}          # 1D array for the differential cross section w.r.t. the kinetic energy of the photoelectrons
#     λe::Array{Cdouble,1}           # 1D array for the electron attenuation length (w.r.t. K_e)


#     # default ctor (it is not really meaningful, it's more for the sake of having a default constructor)
#     function XPSexp() #
#         new(0.0,0.0,0.0,0.0,                                      # experiement parameters
#             0,Array{Cdouble,1}(undef,0),                          # scanned kinetic energy
#             Array{Cdouble,1}(undef,0),Array{Cdouble,1}(undef,0))  # cross section and attenuation length
#     end

#     # ctor: this is physically relevant (as much as the meaning of the parameters used for the model)
#     function XPSexp(ħν_exp::Cdouble,Fν_exp::Cdouble,α_exp::Cdouble,T_exp::Cdouble,θ_exp::Cdouble,Kes::Array{Cdouble,1},σ_exp::Array{Cdouble,1},λ_exp::Array{Cdouble,1})
#         N_exp = length(Kes);
#         if (length(σ_exp)!=N_exp) | (length(λ_exp)!=N_exp)
#             throw("XPSexp: cross section or/and EAL arrays are not the right length")
#         end
#         new(ħν_exp,Fν_exp,α_exp*T_exp,θ_exp,
#             N_exp,Kes,
#             σ_exp,λ_exp)
#     end

#     # cptor
#     function XPSexp(ws::XPSexp) #
#         new(ws.ħν,ws.Fν,ws.αT,ws.θ,ws.Nke,copy(ws.Ke), copy(ws.σ_β),copy(ws.λe))
#     end
# end
