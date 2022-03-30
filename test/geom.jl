## load the packages used in the estimation
# plotting
using PyPlot
fm = PyPlot.matplotlib.font_manager.json_load("/home/mattoz/.cache/matplotlib/fontlist-v310.json")
fm.findfont("serif", rebuild_if_missing=false)
fm.findfont("serif", fontext="afm", rebuild_if_missing=false)
rc("font",family="serif",serif="Computer Modern Roman")
rc("text", usetex=true)
using myPlot

# data manipulation (loading, writing, etc)
using CSV
using Dates
using DataFrames # for loading and saving the results, e.g. CSV.write("measured_psd.csv",DataFrame(PSD_meas'); writeheader=false)
using Printf

# scientific package from the official Julia repositories
using LinearAlgebra
using Statistics
using DSP
using SpecialMatrices
using Polynomials
using StatsBase

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

r = collect(0.0:0.02:μ0);
θ = collect(range(0.0,2π,length=256)) ; #collect(0:0.1:2π);
y = 0.0
K = 256;
Y = 0.0μ0 .+ collect(range(-0.5L,0.5L,length=K));

function dist_polar(r::Cdouble,θ::Cdouble,y::Cdouble,x0::Cdouble,y0::Cdouble,z0::Cdouble)
    R0 = sqrt(z0^2+x0^2);
    θ0 = atan(x0,z0);
    A = -(R0*r*cos(θ-θ0) - r^2)*sqrt(R0^2 - 2R0*r*cos(θ-θ0) + r^2 + (y0-y)^2)/(R0^2 - 2R0*r*cos(θ-θ0) + r^2)
    B = (1.0 + (y0-y)^2/(R0^2 - 2R0*r*cos(θ-θ0) + r^2))*(μ0^2 - (R0^2*r^2*(sin(θ-θ0))^2)/(R0^2 - 2R0*r*cos(θ-θ0) + r^2))
    A+sqrt(B)
end

function dist_polar(r::Array{Cdouble,1},θ::Array{Cdouble,1},y::Cdouble,x0::Cdouble,y0::Cdouble,z0::Cdouble)
    R0 = sqrt(z0^2+x0^2);
    θ0 = atan(x0,z0);
    A = -(R0*r*cos.(θ'.-θ0) .- r.^2).*sqrt.(R0^2 .- 2R0*r*cos.(θ'.-θ0) .+ r.^2 .+ (y0.-y).^2)./(R0^2 .- 2R0*r*cos.(θ'.-θ0) .+ r.^2)
    B = (1.0 .+ (y0-y).^2 ./(R0^2 .- 2R0*r*cos.(θ'.-θ0) .+ r.^2)).*(μ0^2 .- (R0^2*r.^2*(sin.(θ'.-θ0)).^2)./(R0^2 .- 2R0*r*cos.(θ'.-θ0) .+ r.^2))
    A[r.>μ0,:] .= 0.0
    B[r.>μ0,:] .= 0.0
    C = A+sqrt.(B);
    C[C.<0.0] .= 0.0;
    C
end

function dist_polar_simple(r::Array{Cdouble,1},θ::Array{Cdouble,1},y::Cdouble,x0::Cdouble,y0::Cdouble,z0::Cdouble)
    θ0 = atan(x0,z0)
    A = -r*cos.(θ'.-θ0)
    B = (μ0^2 .- (r.^2*(sin.(θ'.-θ0)).^2))
    A[r.>μ0,:] .= 0.0
    B[r.>μ0,:] .= 0.0
    C = A+sqrt.(B);
    C[C.<0.0] .= 0.0;
    C
end

function d_plane_P(x::Cdouble,y::Cdouble,z::Array{Cdouble,1},x0::Cdouble,y0::Cdouble,z0::Cdouble)
    C = -z.*sqrt.(1.0.+((x0.-x).^2 .+(y0.-y).^2)./((z0.-z).^2))
    C[C.<=0.0] .= 0.0;
    C
end


# atan(sqrt(2.0),1.0)
# near the analyzer
x0_near = 1.1*sqrt(2.0);
y0_near = 0.0;
z0_near = 1.1;

DD = dist_polar(r,θ,y,x0_near,y0_near,z0_near);
DD_simple = dist_polar_simple(r,θ,y,x0_near,y0_near,z0_near);

# far away from the analyzer
x0_far = 100.0*sqrt(2.0);
y0_far = 0.0;
z0_far = 100.0;

DD_away = dist_polar(r,θ,y,x0_far,y0_far,z0_far);
DD_simple_away = dist_polar_simple(r,θ,y,x0_far,y0_far,z0_far);


## plot
fig = figure(figsize=[9,5])
# ax = fig.add_axes([0.1,0.1,0.8,0.8],polar=true)
ax1 = subplot(121,polar=true)
pcm1 = ax1.pcolormesh(θ,r,DD,edgecolors="face")
# rc("ytick",color="white")
cax1 = fig.add_axes([0.08, .35, 0.02, 0.3])
cb1 = fig.colorbar(pcm1, orientation="vertical", cax=cax1, shrink=0.6)
cb1.set_label("distance [\$\\mu\$m]", fontsize=10) # , color="white"
# cb1.ax.yaxis.set_tick_params(color="white")
# cb1.outline.set_edgecolor("white")
# rc("ytick",color="black")


ax2 = subplot(122,polar=true)
pcm2 = ax2.pcolormesh(θ,r,DD_simple,edgecolors="face")
# rc("ytick",color="white")
cax2 = fig.add_axes([0.58, .35, 0.02, 0.3])
cb2 = fig.colorbar(pcm2, orientation="vertical", cax=cax2, shrink=0.6)
cb2.set_label("distance [\$\\mu\$m]", fontsize=10) # , color="white"
# cb2.ax.yaxis.set_tick_params(color="white")
# cb2.outline.set_edgecolor("white")
# rc("ytick",color="black")

tight_layout(pad=1.0, w_pad=2.0, h_pad=0.2)

ax1.annotate("a)", xy=(3, 1),  xycoords="data", xytext=(0.0, 1.0), textcoords="axes fraction", color="black",fontsize=14)
ax2.annotate("b)", xy=(3, 1),  xycoords="data", xytext=(0.0, 1.0), textcoords="axes fraction", color="black",fontsize=14)

# savefig("distance_cylinder_near.png")
# savefig("distance_cylinder_near.pdf")



fig = figure(figsize=[9,5])
ax1 = subplot(121,polar=true)
pcm1 = ax1.pcolormesh(θ,r,DD_away,edgecolors="face")
cax1 = fig.add_axes([0.08, .35, 0.02, 0.3])
cb1 = fig.colorbar(pcm1, orientation="vertical", cax=cax1, shrink=0.6)
cb1.set_label("distance [\$\\mu\$m]", fontsize=10)


ax2 = subplot(122,polar=true)
pcm2 = ax2.pcolormesh(θ,r,DD_simple_away,edgecolors="face")
cax2 = fig.add_axes([0.58, .35, 0.02, 0.3])
cb2 = fig.colorbar(pcm2, orientation="vertical", cax=cax2, shrink=0.6)
cb2.set_label("distance [\$\\mu\$m]", fontsize=10)
tight_layout(pad=1.0, w_pad=2.0, h_pad=0.2)

ax1.annotate("a)", xy=(3, 1),  xycoords="data", xytext=(0.0, 1.0), textcoords="axes fraction", color="black",fontsize=14)
ax2.annotate("b)", xy=(3, 1),  xycoords="data", xytext=(0.0, 1.0), textcoords="axes fraction", color="black",fontsize=14)

# savefig("distance_cylinder_far.png")
# savefig("distance_cylinder_far.pdf")



## measurement model
N = 51;
r_surf = collect(range(μ0-5λe,μ0,length=N));

DD_meas_far = dist_polar(r_surf,θ,y,x0_far,y0_far,z0_far);
# DD_meas = dist_polar_simple(r_surf,θ,y,x0_far,y0_far,z0_far);
DD_meas_near = dist_polar(r_surf,θ,y,x0_near,y0_near,z0_near);

intDD_far = r_surf.*exp.(-DD_meas_far./λe);
intDD_far_θ = sum(intDD_far,dims=2);
intDD_near = r_surf.*exp.(-DD_meas_near./λe);
intDD_near_θ = sum(intDD_near,dims=2);


## compare to simple exp model (planar case)
relative_eal_far  = log.(exp.(-(μ0.-r_surf)/λe))./log.((intDD_far_θ/intDD_far_θ[end]));
relative_eal_near = log.(exp.(-(μ0.-r_surf)/λe))./log.((intDD_near_θ/intDD_near_θ[end]));
r_eal_min_far,r_eal_max_far = extrema(relative_eal_far[1:end-1])
r_eal_min_near,r_eal_max_near = extrema(relative_eal_near[1:end-1])

figure();
plot(r_surf,relative_eal_far, label="far")
plot(r_surf,relative_eal_near, label="near")
legend()
xlabel("radius [\$\\mu m\$]")
ylabel("log ratio: \$\\frac{\\log(\\frac{I^{\\mathrm{planar}}}{I^{\\mathrm{planar}}_0} )}{\\log(\\frac{I^{\\mathrm{cylindrical}}}{I^{\\mathrm{cylindrical}}_0})}\$")

## show variation in aparent eal using planar model to approximate the cylindrical model
figure();
plot(r_surf,intDD_far_θ/intDD_far_θ[end], label="cylinder: far \$\\lambda_e\$", color="tab:blue");
plot(r_surf,intDD_near_θ/intDD_near_θ[end], label="cylinder: near \$\\lambda_e\$", color="tab:green");
plot(r_surf,exp.(-(μ0.-r_surf)/λe), label="planar approximation \$\\lambda_e\$", color="tab:orange");
s = @sprintf "planar approximation limits (far) \$[%.2f\\lambda_e,%.2f\\lambda_e]\$" r_eal_min_far r_eal_max_far
fill_between(r_surf,exp.(-(μ0.-r_surf)/(r_eal_min_far*λe)),exp.(-(μ0.-r_surf)/(r_eal_max_far*λe)),alpha=0.5,color="tab:red",label=s)
s = @sprintf "planar approximation limits (near) \$[%.2f\\lambda_e,%.2f\\lambda_e]\$" r_eal_min_near r_eal_max_near
fill_between(r_surf,exp.(-(μ0.-r_surf)/(r_eal_min_near*λe)),exp.(-(μ0.-r_surf)/(r_eal_max_near*λe)),alpha=0.5,color="tab:pink",label=s)
legend()
xlabel("radius [\$\\mu m\$]")
ylabel("normalized gain")
# savefig("plans_vs_cylinder.png")
# savefig("plans_vs_cylinder.pdf")


## compare the cylindrical acquisition gain density in the near and far cases
fig = figure(figsize=[9,5])
# near analyzer
ax1 = subplot(121,polar=true)
ax1.set_rlabel_position(-22.5)
ax1.set_rticks([0.97, 0.98, 0.99, 1.0])
ax1.plot(atan(x0_far,z0_far)*ones(Cdouble,2),[0.97; 1.0], color="red",label="\$\\theta\\simeq54.7\$")
pcm1 = ax1.pcolormesh(θ,r_surf,intDD_near,edgecolors="face")
ylim(0.97,1.0)
ax1.legend(loc="lower left", bbox_to_anchor=(.5 + cos(atan(x0_far,z0_far))/2, .5 + sin(atan(x0_far,z0_far))/2))
cax1 = fig.add_axes([0.15, .37, 0.02, 0.25])
cb1 = fig.colorbar(pcm1, orientation="vertical", cax=cax1, shrink=0.6)
cb1.set_label("gain [a.u.]", fontsize=10) # , color="white"

#far analyzer
ax2 = subplot(122,polar=true)
ax2.set_rlabel_position(-22.5)
ax2.set_rticks([0.97, 0.98, 0.99, 1.0])
ax2.plot(atan(x0_far,z0_far)*ones(Cdouble,2),[0.97; 1.0], color="red",label="\$\\theta\\simeq54.7\$")
pcm2 = ax2.pcolormesh(θ,r_surf,intDD_far,edgecolors="face")
ylim(0.97,1.0)
ax2.legend(loc="lower left", bbox_to_anchor=(.5 + cos(atan(x0_far,z0_far))/2, .5 + sin(atan(x0_far,z0_far))/2))
cax2 = fig.add_axes([0.65, .37, 0.02, 0.25])
cb2 = fig.colorbar(pcm2, orientation="vertical", cax=cax2, shrink=0.6)
cb2.set_label("gain [a.u.]", fontsize=10) # , color="white"

# annotation and squeezing
tight_layout(pad=1.0, w_pad=2.0, h_pad=0.2)
ax1.annotate("a)", xy=(3, 1),  xycoords="data", xytext=(0.0, 1.0), textcoords="axes fraction", color="black",fontsize=14)
ax2.annotate("b)", xy=(3, 1),  xycoords="data", xytext=(0.0, 1.0), textcoords="axes fraction", color="black",fontsize=14)

# savefig("gain_cylinder_near.png")
# savefig("gain_cylinder_near.pdf")


## compute the complete model

J = 257;
θ0_far  = atan(x0_far,z0_far);
θ_far   = collect(range(θ0_far-π/2.0,θ0_far+π/2.0,length=J));
θ0_near = atan(x0_near,z0_near);
θ_near  = collect(range(θ0_near-π/2.0,θ0_near+π/2.0,length=J));

function cylinder_gain_H(r::Array{Cdouble,1},θ::Array{Cdouble,1},y::Array{Cdouble,1},x0::Cdouble,y0::Cdouble,z0::Cdouble,μ0::Cdouble,λe::Cdouble)
    # compute the elementary integrals
    Arn = 0.5*r.*[r[2]-r[1]; r[3:end]-r[1:end-2]; r[end]-r[end-1]]; # rdr
    Aθj = 0.5*[θ[2]-θ[1]; θ[3:end]-θ[1:end-2]; θ[end]-θ[end-1]];    # dθ
    Ayk = 0.5*[y[2]-y[1]; y[3:end]-y[1:end-2]; y[end]-Y[end-1]];    # dy

    #compute the model
    H_rθy = zeros(Cdouble,length(r),length(θ),length(y));
    for k in 1:length(y)
        H_rθy[:,:,k] = exp.(-d_cylinder_P(r,θ,y[k],x0,y0,z0,μ0)/λe)
    end

    H_r  = zeros(Cdouble,length(r));
    for n in 1:length(r)
        H_r[n] = Arn[n]*Aθj'*H_rθy[n,:,:]*Ayk
    end
    H_r,H_rθy,Arn,Aθj,Ayk
end

#TODO: check if the results are OK
H_r,H_rθy,Arn,Aθj,Ayk = cylinder_gain_H(r_surf,θ_far,Y,x0_far,y0_far,z0_far,μ0,λe)
H_r,H_rθy,Arn,Aθj,Ayk = cylinder_gain_H(r_surf,θ_near,Y,x0_near,y0_near,z0_near,μ0,λe)

Arn = 0.5*[r_surf[2]-r_surf[1]; r_surf[3:end]-r_surf[1:end-2]; r_surf[end]-r_surf[end-1]];
Aθj_far = 0.5*[θ_far[2]-θ_far[1]; θ_far[3:end]-θ_far[1:end-2]; θ_far[end]-θ_far[end-1]];
Aθj_near = 0.5*[θ_near[2]-θ_near[1]; θ_near[3:end]-θ_near[1:end-2]; θ_near[end]-θ_near[end-1]];
Ayk = 0.5*[Y[2]-Y[1]; Y[3:end]-Y[1:end-2]; Y[end]-Y[end-1]];

H_far = zeros(Cdouble,N,J,K);
H_near = zeros(Cdouble,N,J,K);
for k in 1:K
    H_far[:,:,k] = exp.(-dist_polar(r_surf,θ_far,Y[k],x0_far,y0_far,z0_far)/λe)
    H_near[:,:,k] = exp.(-dist_polar(r_surf,θ_near,Y[k],x0_near,y0_near,z0_near)/λe)
end
H_far = H_far.*r_surf;
H_near = H_near.*r_surf;

H_n_far  = zeros(Cdouble,N);
H_n_near = zeros(Cdouble,N);
for n in 1:N
    H_n_far[n] = Arn[n]*Aθj_far'*H_far[n,:,:]*Ayk
    H_n_near[n] = Arn[n]*Aθj_near'*H_near[n,:,:]*Ayk
end


# compare to planar model (need to normalize)
H_z = Arn.*exp.(-(μ0.-r_surf)/λe);
H_z_far_min = Arn.*exp.(-(μ0.-r_surf)/(0.5436*λe));
H_z_far_max = Arn.*exp.(-(μ0.-r_surf)/(0.9193*λe));
H_z_near_min = Arn.*exp.(-(μ0.-r_surf)/(0.4034*λe));
H_z_near_max = Arn.*exp.(-(μ0.-r_surf)/(0.8364*λe));


figure();
plot(r_surf,H_n_far/maximum(H_n_far), label="cylinder: far \$\\lambda_e\$", color="tab:blue");
plot(r_surf,H_n_near/maximum(H_n_near), label="cylinder: near \$\\lambda_e\$", color="tab:green");
plot(r_surf,H_z/maximum(H_z), label="planar approximation \$\\lambda_e\$", color="tab:orange");
s = @sprintf "planar approximation limits (far)" #  \$[%.2f\\lambda_e,%.2f\\lambda_e]\$" r_eal_min_far r_eal_max_far
fill_between(r_surf,H_z_far_min/maximum(H_z_far_min),H_z_far_max/maximum(H_z_far_max),alpha=0.5,color="tab:red",label=s)
s = @sprintf "planar approximation limits (near)" #  \$[%.2f\\lambda_e,%.2f\\lambda_e]\$" r_eal_min_near r_eal_max_near
fill_between(r_surf,H_z_near_min/maximum(H_z_near_min),H_z_near_max/maximum(H_z_near_max),alpha=0.5,color="tab:pink",label=s)
legend()
xlabel("radius [\$\\mu m\$]")
ylabel("normalized gain")



# simulate some data (one point in the kinetic energy spectrum for four different concentration profiles)
ρA_1 = logistic.(1000.0reverse(μ0.-r_surf).-2.0,0.0,1.0,2.0);
ρA_2 = logistic.(1000.0reverse(μ0.-r_surf).-2.0,0.0,1.0,2.0) .+ 2.0exp.(-(1000.0reverse(μ0.-r_surf).-1.0).^2. /(2.0*0.25^2));
ρA_3 = logistic.(1000.0reverse(μ0.-r_surf).-2.0,0.0,1.0,2.0) .+ exp.(-(1000.0reverse(μ0.-r_surf).-1.5).^2. /(2.0*0.5^2));
ρA_4 = exp.(-(1000.0reverse(μ0.-r_surf).-2.5).^2. /(2.0*0.5^2));

# figure();
# plot(r_surf,reverse(ρA_1))
# plot(r_surf,reverse(ρA_2))
# plot(r_surf,reverse(ρA_3))
# plot(r_surf,reverse(ρA_4))

# for each cases (cylindrical near and far, and planar), simulate the acquisition of one datum for each profile
M_far  = [H_n_far'*reverse(ρA_1); H_n_far'*reverse(ρA_2); H_n_far'*reverse(ρA_3); H_n_far'*reverse(ρA_4)];
M_near = [H_n_near'*reverse(ρA_1); H_n_near'*reverse(ρA_2); H_n_near'*reverse(ρA_3); H_n_near'*reverse(ρA_4)];
M_z    = [H_z'*reverse(ρA_1); H_z'*reverse(ρA_2); H_z'*reverse(ρA_3); H_z'*reverse(ρA_4)];

figure();
scatter([1; 2; 3; 4],M_far)
scatter([1; 2; 3; 4],M_near)
scatter([1; 2; 3; 4],M_z)

figure(figsize=[10,5])
ax1 = subplot(121)
plot(r_surf,H_n_far/maximum(H_n_far), label="cylinder: far \$\\lambda_e\$", color="tab:blue");
plot(r_surf,H_n_near/maximum(H_n_near), label="cylinder: near \$\\lambda_e\$", color="tab:green");
plot(r_surf,H_z/maximum(H_z), label="planar approximation \$\\lambda_e\$", color="tab:orange");
s = @sprintf "planar approximation limits (far)" #  \$[%.2f\\lambda_e,%.2f\\lambda_e]\$" r_eal_min_far r_eal_max_far
fill_between(r_surf,H_z_far_min/maximum(H_z_far_min),H_z_far_max/maximum(H_z_far_max),alpha=0.5,color="tab:red",label=s)
s = @sprintf "planar approximation limits (near)" #  \$[%.2f\\lambda_e,%.2f\\lambda_e]\$" r_eal_min_near r_eal_max_near
fill_between(r_surf,H_z_near_min/maximum(H_z_near_min),H_z_near_max/maximum(H_z_near_max),alpha=0.5,color="tab:pink",label=s)
legend()
xlabel("radius [\$\\mu m\$]")
ylabel("normalized gain")

ax2 = subplot(122)
scatter(["\$\\rho_1\$"; "\$\\rho_2\$"; "\$\\rho_3\$"; "\$\\rho_4\$"],M_far./M_z, label="cylinder far vs planar", color="tab:blue")
scatter(["\$\\rho_1\$"; "\$\\rho_2\$"; "\$\\rho_3\$"; "\$\\rho_4\$"],M_near./M_z, label="cylinder near vs planar", color="tab:green")
xlabel("concentration profiles")
ylabel("relative acquisition: \$\\frac{I_{\\mathrm{cylinder}}}{I_{\\mathrm{plane}}}\$")
ylim(0.0)
legend()

tight_layout(pad=1.0, w_pad=2.0, h_pad=0.2)
ax2.annotate("a)", xy=(3, 1),  xycoords="data", xytext=(-1.27, 0.99), textcoords="axes fraction", color="black",fontsize=14)
ax2.annotate("b)", xy=(3, 1),  xycoords="data", xytext=(-0.07, 0.99), textcoords="axes fraction", color="black",fontsize=14)

# savefig("plans_vs_cylinder.png")
# savefig("plans_vs_cylinder.pdf")

figure()
scatter([1; 2; 3; 4],(M_far./M_z)./(M_near./M_z))


##
## planar model with extent over some area, not just one point in space
##

Nx = 100;
Ny = K;
Nz = N;
x_far   = collect(range(-μ0,μ0,length=Nx));
# y_far   = collect(range(y0_far-L/2.0,y0_far+L/2.0,length=Ny));
x_near  = collect(range(-μ0,μ0,length=Nx));
# y_near  = collect(range(y0_near-L/2.0,y0_near+L/2.0,length=Ny));

z_surf = collect(range(-5λe,0.0,length=N));

Azn      = 0.5*[z_surf[2]-z_surf[1]; z_surf[3:end]-z_surf[1:end-2]; z_surf[end]-z_surf[end-1]];
Axj_far  = 0.5*[x_far[2]-x_far[1]; x_far[3:end]-x_far[1:end-2]; x_far[end]-x_far[end-1]];
Axj_near = 0.5*[x_near[2]-x_near[1]; x_near[3:end]-x_near[1:end-2]; x_near[end]-x_near[end-1]];
Ayk = 0.5*[Y[2]-Y[1]; Y[3:end]-Y[1:end-2]; Y[end]-Y[end-1]];

H_z_far = zeros(Cdouble,Nz,Nx,Ny);
H_z_near = zeros(Cdouble,Nz,Nx,Ny);

for j in 1:Nx
    for k in 1:Ny
        H_z_far[:,j,k] = exp.(-d_plane_P(x_far[j],Y[k],z_surf,x0_far,y0_far,z0_far)/λe)
        H_z_near[:,j,k] = exp.(-d_plane_P(x_near[j],Y[k],z_surf,x0_near,y0_near,z0_near)/λe)
    end
end

H_z_far_int  = zeros(Cdouble,Nz);
H_z_near_int = zeros(Cdouble,Nz);

for n in 1:Nz
    H_z_far_int[n]  = Axj_far'*H_z_far[n,:,:]*Ayk
    H_z_near_int[n] = Axj_near'*H_z_near[n,:,:]*Ayk
end

H_z_far_int = H_z_far_int.*Azn
H_z_near_int = H_z_near_int.*Azn
H_z_z = Arn.*exp.(z_surf/λe);
H_z_z_0 = Arn.*exp.(z_surf/(0.58λe));

figure()
plot(z_surf,H_z_far_int/maximum(H_z_far_int))
plot(z_surf,H_z_near_int/maximum(H_z_near_int))
plot(z_surf,H_z_z/maximum(H_z_z))
plot(z_surf,H_z_z_0/maximum(H_z_z_0))


function plane_gain_H(x::Array{Cdouble,1},y::Array{Cdouble,1},z::Array{Cdouble,1},x0::Cdouble,y0::Cdouble,z0::Cdouble,λe::Cdouble)
    # compute the elementary integrals
    Azn = 0.5*[z[2]-z[1]; z[3:end]-z[1:end-2]; z[end]-z[end-1]]; # dz
    Axj = 0.5*[x[2]-x[1]; x[3:end]-x[1:end-2]; x[end]-x[end-1]]; # dx
    Ayk = 0.5*[y[2]-y[1]; y[3:end]-y[1:end-2]; y[end]-Y[end-1]]; # dy

    #compute the model
    H_zxy = zeros(Cdouble,length(z),length(x),length(y));

    for j in 1:length(x)
        for k in 1:length(y)
            H_zxy[:,j,k] = exp.(-d_plane_P(x[j],y[k],z,x0,y0,z0)/λe)
        end
    end

    H_z  = zeros(Cdouble,length(z));
    for n in 1:length(z)
        H_z[n]  = Azn[n]*Axj'*H_z_far[n,:,:]*Ayk
    end

    H_z,H_zxy,Azn,Axj,Ayk
end
