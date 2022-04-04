## load the packages used in the estimation
# plotting
using PyPlot
fm = PyPlot.matplotlib.font_manager.json_load("/home/matthew/.cache/matplotlib/fontlist-v310.json") # TODO: look for the path automatically
# fm = PyPlot.matplotlib.font_manager.json_load("/home/mattoz/.cache/matplotlib/fontlist-v310.json")
# fm.findfont("serif", rebuild_if_missing=false)
# fm.findfont("serif", fontext="afm", rebuild_if_missing=false)
rc("font",family="serif",serif="Computer Modern Roman")
rc("text", usetex=true)
using myPlot

# data manipulation (loading, writing, etc)
using Printf

# scientific package from the official Julia repositories
# using LinearAlgebra
# using Statistics
# using DSP
# using SpecialMatrices
# using Polynomials
# using StatsBase

# implemented scientific packages
using utilsFun  # for the softMax functions

# modeling XPS
using XPSpack
using XPSinv

# Dict{Int64,XPSsetup}();

# the geometry structures assumes that
#  - the photon beam is along the z axis
#  - the liquid microjet is along the y axis (or the droplet axis motion, or, if no clear axis arise naturally, any direction orthogonal to the z axis)
#  - the x axis is the remaining axis that makes Oxyz an orthogonal frame of reference
#  - the reference of axis, namely O, is taken either at the center or on the surface of the sample
# maybe but not sure: the model derived from the geomStructs assumes that the incident light is orthogonal to the surface of the sample

# # finger struct
# mutable struct fingerGeom
#     # coordinates of the analyzer's apperture (P)
#     x0::Cdouble
#     y0::Cdouble
#     z0::Cdouble

#     # coordinate of the finger at the sample's surface
#     x::Cdouble
#     y::Cdouble

#     # probing depths along the finger (since the sampling can be non-uniform, I decide to keep an array instead of a range)
#     z::Array{Cdouble,1}

#     function fingerGeom()
#         new(0.0,0.0,0.0, 0.0,0.0,Array{Cdouble,1}(undef,0))
#     end
#     function fingerGeom(x0_::Cdouble,y0_::Cdouble,z0_::Cdouble,x_::Cdouble,y_::Cdouble,z_::Array{Cdouble,1})
#         new(x0_,y0_,z0_,x_,y_,z_)
#     end
#     function fingerGeom(ws::fingerGeom)
#         new(ws.x0,ws.y0,ws.z0,ws.x,ws.y,ws.z)
#     end
# end

# # plane struct: assumes a rectangular cuboid discretization volume
# mutable struct planeGeom
#     # coordinates of the analyzer's apperture (P)
#     x0::Cdouble
#     y0::Cdouble
#     z0::Cdouble

#     # coordinates of the illuminated area at the surface of the sample (the covered area may be bigger or smaller than the sample's dimension: the derived model should deal with that with the beam-light's profile)
#     x::Array{Cdouble,1}
#     y::Array{Cdouble,1}

#     # probing depths
#     z::Array{Cdouble,1}

#     function planeGeom()
#         new(0.0,0.0,0.0, Array{Cdouble,1}(undef,0),Array{Cdouble,1}(undef,0),Array{Cdouble,1}(undef,0))
#     end
#     function planeGeom(x0_::Cdouble,y0_::Cdouble,z0_::Cdouble,x_::Array{Cdouble,1},y_::Array{Cdouble,1},z_::Array{Cdouble,1})
#         new(x0_,y0_,z0_,x_,y_,z_)
#     end
#     function planeGeom(ws::planeGeom)
#         new(ws.x0,ws.y0,ws.z0,ws.x,ws.y,ws.z)
#     end
# end

# # cylinder struct
# mutable struct cylinderGeom
#     # coordinates of the analyzer's apperture (P)
#     x0::Cdouble
#     y0::Cdouble
#     z0::Cdouble

#     # radius of the sample
#     μ0::Cdouble 

#     # probing depths (distance from the symmetry axis of the sample)
#     r::Array{Cdouble,1}

#     # coordinates of the illuminated area at the surface of the sample in polar coordinates (the covered area may be bigger or smaller than the sample's dimension: the derived model should deal with that with the beam-light's profile)
#     θ::Array{Cdouble,1}
#     y::Array{Cdouble,1}

#     function cylinderGeom()
#         new(0.0,0.0,0.0,0.0, Array{Cdouble,1}(undef,0),Array{Cdouble,1}(undef,0),Array{Cdouble,1}(undef,0))
#     end
#     function cylinderGeom(x0_::Cdouble,y0_::Cdouble,z0_::Cdouble,μ0_::Cdouble,r_::Array{Cdouble,1},θ_::Array{Cdouble,1},y_::Array{Cdouble,1})
#         new(x0_,y0_,z0_,μ0_,r_,θ_,y_)
#     end
#     function cylinderGeom(ws::cylinderGeom)
#         new(ws.x0,ws.y0,ws.z0,ws.μ0,ws.r,ws.θ,ws.y)
#     end
# end


# # acquisition struct
# mutable struct XPSacq
#     # main setup number
#     ħν::Cdouble                    # photon energy (array of length Nν)
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


## acquisition setup
ħν = 900.0;
μKe = ħν-285.0;
α = 1.0
T = 1.0
Fν = 1.0;
Nke = 200;
Ke = collect(range(μKe-2.0,μKe+2.0,length=Nke));
dKe = Ke[2]-Ke[1]
σν0 = 0.6;
σν = σν0*((0.7/sqrt(2π*0.2^2))*exp.(-0.5*(Ke.-(μKe-0.5)).^2/0.2^2) .+ (0.3/sqrt(2π*0.5^2))*exp.(-0.5*(Ke.-(μKe+0.5)).^2/0.5^2));
λe0 = 0.002;

wsAcq = XPSacq(ħν,μKe,α,T,Fν,Ke,σν,λe0);

## geometry setup
k0 = 5;
Nr = 51;
Nθ = 256;
Ny = 256;
μ0 = 100.0;
L = 200.0;
x0 = sqrt(2.0)*100.0#μ0
y0 = 0.0;
z0 = 100.0#μ0
r = collect(range(μ0-k0*λe0,μ0,length=Nr))
θ0 = atan(x0,z0)
θ = collect(range(θ0-π/2.0,θ0+π/2.0,Nθ));
y = collect(range(-L/2.0,L/2.0,length=Ny));

wsGeom = cylinderGeom(x0,y0,z0,μ0,r,θ,y)

# TODO: create a model from the given elements
Hr,Hrθy,Arn,Aθj,Ayk = cylinder_gain_H(r,θ,y,x0,y0,z0,μ0,λe0);



fig,ax,pcm,cax,cb = imshowDataPolar(1,r,θ,Hrθy[:,:,128];cb_ax_loc=(0.25, .37));
# ax.set_rticks([99.97, 99.98, 99.99, 100.0])
ax.set_ylim(99.97,100.0)



function acqModel(wsAcq::XPSacq,wsGeom::cylinderGeom)
    Hr,Hrθy,_,_,_ = cylinder_gain_H(wsGeom.r,wsGeom.θ,wsGeom.y,wsGeom.x0,wsGeom.y0,wsGeom.z0,wsGeom.μ0,wsAcq.λe);
    wsAcq.T*wsAcq.α*wsAcq.Fν*σν*Hr'
end
