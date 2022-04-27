## load the packages used in the estimation
# plotting
using PyPlot
fm = PyPlot.matplotlib.font_manager.json_load("/home/matthew/.cache/matplotlib/fontlist-v310.json")
# fm.findfont("serif", rebuild_if_missing=false)
# fm.findfont("serif", fontext="afm", rebuild_if_missing=false)
rc("font",family="serif",serif="Computer Modern Roman")
rc("text", usetex=true)

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



##
## distances
##

μ0 = 1.0; # radius of the microjet
λe = 2.0e-3μ0; # EAL
L = 3μ0;

N = 51;
K = 256;
J = 257
r = collect(range(0.0,μ0,length=N));
θ = collect(range(0.0,2π,length=J)) ; #collect(0:0.1:2π);
y = 0.0
Y = 0.0μ0 .+ collect(range(-0.5L,0.5L,length=K));


# magic angle: atan(sqrt(2.0),1.0)
# near the analyzer
x0_near = 1.1*sqrt(2.0);
y0_near = 0.0;
z0_near = 1.1;

# far away from the analyzer
x0_far = 100.0*sqrt(2.0);
y0_far = 0.0;
z0_far = 100.0;

## compute and plot distances for two cases (near and far), and compare the full model and its approximation
include("distance.jl")


##
## cylindrical model and finger model
##

## model specific discretization
r_surf = collect(range(μ0-5λe,μ0,length=N));
θ0_far  = atan(x0_far,z0_far);
θ_far   = collect(range(θ0_far-π/2.0,θ0_far+π/2.0,length=J));
θ0_near = atan(x0_near,z0_near);
θ_near  = collect(range(θ0_near-π/2.0,θ0_near+π/2.0,length=J));

## compute and plot the models
include("finger_and_cylinder_model.jl")



##
## planar model with extent over some area, not just one point in space
##

Nx = 100;
Ny = K;
Nz = N;
x_far   = collect(range(-μ0,μ0,length=Nx));
x_near  = collect(range(-μ0,μ0,length=Nx));
z_surf = collect(range(-5λe,0.0,length=N));

## compute models and plot
include("finger_and_planar_model.jl")



##
## acquisition simulation
##

# simulate some data (one point in the kinetic energy spectrum for four different concentration profiles)
ρA_1 = logistic.(1000.0reverse(μ0.-r_surf).-2.0,0.0,1.0,2.0);
ρA_2 = logistic.(1000.0reverse(μ0.-r_surf).-2.0,0.0,1.0,2.0) .+ 2.0exp.(-(1000.0reverse(μ0.-r_surf).-1.0).^2. /(2.0*0.25^2));
ρA_3 = logistic.(1000.0reverse(μ0.-r_surf).-2.0,0.0,1.0,2.0) .+ exp.(-(1000.0reverse(μ0.-r_surf).-1.5).^2. /(2.0*0.5^2));
ρA_4 = exp.(-(1000.0reverse(μ0.-r_surf).-2.5).^2. /(2.0*0.5^2));

include("geom_acquisition.jl")
