## load the packages used in the estimation
# plotting
using PyPlot
rc("text", usetex=true)

# data manipulation (loading, writing, etc)
using Printf

# implemented scientific packages
using utilsFun  # for the softMax functions

# modeling XPS
using XPSpack


##
## distances
##

μ0 = 10.0; # radius of the microjet
δr = 2.0e-3 # transition to vacuum layer thickness (let's set about 1 nm)
λe = 2.0e-3; # μ0; # EAL
L = 3μ0;

N = 51;
K = 256;
J = 257
r = collect(range(0.0,μ0+δr,length=N));
θ = collect(range(0.0,2π,length=J)) ; #collect(0:0.1:2π);
y = 0.0
Y = 0.0μ0 .+ collect(range(-0.5L,0.5L,length=K));


# magic angle: atan(sqrt(2.0),1.0)
# near the analyzer
x0_near = 21.0*sqrt(2.0);
y0_near = 0.0;
z0_near = 21.0;

# far away from the analyzer
x0_far = 200.0*sqrt(2.0);
y0_far = 0.0;
z0_far = 200.0;

## compute and plot distances for two cases (near and far), and compare the full model and its approximation
include("distance_vapor.jl")


##
## cylindrical model and finger model
##

## model specific discretization
r_surf = collect(range(μ0-5λe,μ0+δr,length=N));
θ0_far  = atan(x0_far,z0_far);
θ_far   = collect(range(θ0_far-π/2.0,θ0_far+π/2.0,length=J));
θ0_near = atan(x0_near,z0_near);
θ_near  = collect(range(θ0_near-π/2.0,θ0_near+π/2.0,length=J));

## compute and plot the models
include("finger_and_cylinder_model_vapor.jl")



##
## acquisition simulation
##

# simulate some data (one point in the kinetic energy spectrum for four different concentration profiles)
ρ0 = 1.0
ρ_vac = 0.0
r_th  = 2.346; # 2.0 #
ρA_1 = logistic.(1000.0reverse(δr.+μ0.-r_surf).-r_th,ρ_vac,ρ0,2.0);
# ρA_2 = logistic.(1000.0reverse(δr.+μ0.-r_surf).-r_th,ρ_vac,ρ0,2.0) .+ 2.0exp.(-(1000.0reverse(μ0.-r_surf).-1.0).^2. /(2.0*0.25^2));
ρA_2 = logistic.(1000.0reverse(δr.+μ0.-r_surf).-r_th,ρ_vac,ρ0,2.0) .+ 2.0exp.(-(1000.0reverse(μ0.-r_surf).-0.0).^2. /(2.0*0.25^2));
# ρA_3 = logistic.(1000.0reverse(δr.+μ0.-r_surf).-r_th,ρ_vac,ρ0,2.0) .+ exp.(-(1000.0reverse(μ0.-r_surf).-1.5).^2. /(2.0*0.5^2));
ρA_3 = logistic.(1000.0reverse(δr.+μ0.-r_surf).-r_th,ρ_vac,ρ0,2.0) .+ exp.(-(1000.0reverse(μ0.-r_surf).-0.5).^2. /(2.0*0.5^2));
# ρA_4 = exp.(-(1000.0reverse(δr.+μ0.-r_surf).-2.5).^2. /(2.0*0.5^2));
ρA_4 = exp.(-(1000.0(reverse(δr.+μ0.-r_surf).-δr)).^2. /(2.0*0.5^2));

figure(); 
plot(reverse(μ0.-r_surf),ρA_1)
plot(reverse(μ0.-r_surf),ρA_2)
plot(reverse(μ0.-r_surf),ρA_3)
plot(reverse(μ0.-r_surf),ρA_4)

include("geom_acquisition_vapor.jl")


