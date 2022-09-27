## load the packages used in the estimation
# plotting
using PyPlot
rc("text", usetex=true)

# data manipulation (loading, writing, etc)
using Printf

# modeling XPS
using XPSpack



##
## distances
##

μ0 = 20.0; # radius of the microjet
λe = 2.0e-3; # μ0; # EAL
k0 = 5.0;
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
x0_near = 21.0*sqrt(2.0);
y0_near = 0.0;
z0_near = 21.0;

# far away from the analyzer
x0_far = 200.0*sqrt(2.0);
y0_far = 0.0;
z0_far = 200.0;

## compute and plot distances for two cases (near and far), and compare the full model and its approximation
include("distance.jl")


##
## cylindrical model and finger model
##

## model specific discretization
r_surf = collect(range(μ0-k0*λe,μ0,length=N));
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
z_surf = collect(range(-k0*λe,0.0,length=N));

## compute models and plot
include("finger_and_planar_model.jl")


##
## spherical model
##
# r = collect(range(μ0-k0*λe0,μ0,length=Nr)); # +δr
# θ0 = atan(x0,z0);
# φ0 = acos(y0/sqrt(x0^2 + y0^2 + z0^2));
# φ = collect(range(φ0-π/2.0,φ0+π/2.0,Nφ)); 
# θ = collect(range(θ0-π/2.0,θ0+π/2.0,Nθ));
Nθ = K;
Nφ = J;
φ_sphere = collect(range(0.0,π,Nφ));    # the polar angle (different notations compared with )
θ_sphere = collect(range(0.0,2.0π,Nθ)); # the azimuthal angle

include("finger_and_sphere_model.jl")


##
## acquisition simulation
##

# simulate some data (one point in the kinetic energy spectrum for four different concentration profiles)
ρA_1 = logistic.(1000.0reverse(μ0.-r_surf).-2.0,0.0,1.0,2.0);
ρA_2 = logistic.(1000.0reverse(μ0.-r_surf).-2.0,0.0,1.0,2.0) .+ 2.0exp.(-(1000.0reverse(μ0.-r_surf).-1.0).^2. /(2.0*0.25^2));
ρA_3 = logistic.(1000.0reverse(μ0.-r_surf).-2.0,0.0,1.0,2.0) .+ exp.(-(1000.0reverse(μ0.-r_surf).-1.5).^2. /(2.0*0.5^2));
ρA_4 = exp.(-(1000.0reverse(μ0.-r_surf).-2.5).^2. /(2.0*0.5^2));

include("geom_acquisition.jl")








# ##
# ## plot distance and geometry factor in the away case
# ##


# ## plot
# fig = figure(figsize=[10,5])
# ax1 = subplot(121,polar=true)
# ax1.set_rticks([μ0/4, 2μ0/4, 3μ0/4, μ0])
# yticks(fontsize=12)
# xticks(fontsize=12)
# # ax1.set_rlabel_position(-22.5)
# ax1.set_rlabel_position(234.7)
# ax1.plot(atan(x0_far,z0_far)*ones(Cdouble,2),[0.0; μ0], color="red",label="\$\\theta\\simeq54.7\$")
# pcm1 = ax1.pcolormesh(θ,r,DD_away,edgecolors="face")
# cax1 = fig.add_axes([0.08, .35, 0.02, 0.3])
# cb1 = fig.colorbar(pcm1, orientation="vertical", cax=cax1, shrink=0.6)
# cb1.set_label("distance [\$\\mu\$m]", fontsize=12)
# cb1.ax.tick_params(labelsize=12)
# ax1.legend(loc="lower left", bbox_to_anchor=(.5 + cos(atan(x0_far,z0_far))/2, .5 + sin(atan(x0_far,z0_far))/2),fontsize=12)


# #far analyzer
# ax2 = subplot(122,polar=true)
# # ax2.set_rlabel_position(-22.5)
# ax2.set_rlabel_position(234.7)
# ax2.set_rticks([μ0-0.03, μ0-0.02, μ0-0.01, μ0])
# yticks(fontsize=12)
# xticks(fontsize=12)
# ax2.plot(atan(x0_far,z0_far)*ones(Cdouble,2),[μ0-0.03; μ0], color="red",label="\$\\theta\\simeq54.7\$")
# pcm2 = ax2.pcolormesh(θ_far,r_surf,H_rθy_far[:,:,128],edgecolors="face")
# ylim(μ0-0.03,μ0)
# ax2.legend(loc="lower left", bbox_to_anchor=(.5 + cos(atan(x0_far,z0_far))/2, .5 + sin(atan(x0_far,z0_far))/2),fontsize=12)
# cax2 = fig.add_axes([0.65-0.07, .35, 0.02, 0.3])
# cb2  = fig.colorbar(pcm2, orientation="vertical", cax=cax2, shrink=0.6)
# cb2.set_label("gain [a.u.]", fontsize=12) # , color="white"
# cb2.ax.tick_params(labelsize=12)



# tight_layout(pad=1.0, w_pad=2.0, h_pad=0.2)

# ax1.annotate("a)", xy=(3, 1),  xycoords="data", xytext=(0.0, 1.0), textcoords="axes fraction", color="black",fontsize=14)
# ax1.annotate("b)", xy=(3, 1),  xycoords="data", xytext=(1.3, 1.0), textcoords="axes fraction", color="black",fontsize=14)

# fig.savefig("distance_and_gain_far_20_microns.png")
# fig.savefig("distance_and_gain_far_20_microns.pdf")

