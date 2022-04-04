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
