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

r = collect(0.0:0.02:μ0);
θ = collect(range(0.0,2π,length=256)) ; #collect(0:0.1:2π);
y = 0.0


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

# near the analyzer
x0_near = 1.1;
y0_near = 0.0;
z0_near = 1.1;

DD = dist_polar(r,θ,y,x0_near,y0_near,z0_near);
DD_simple = dist_polar_simple(r,θ,y,x0_near,y0_near,z0_near);

# far away from the analyzer
x0_far = 100.0;
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
r_surf = collect(range(μ0-5λe,μ0,length=51));

# DD_meas = dist_polar(r_surf,θ,y,x0_far,y0_far,z0_far);
# DD_meas = dist_polar_simple(r_surf,θ,y,x0_far,y0_far,z0_far);
DD_meas = dist_polar(r_surf,θ,y,x0_near,y0_near,z0_near);

intDD = r_surf.*exp.(-DD_meas./λe);
intDD_θ = sum(intDD,dims=2);


relative_eal = log.(exp.(-(μ0.-r_surf)/λe))./log.((intDD_θ/intDD_θ[end]))
figure(); plot(r_surf,relative_eal)
r_eal_min,r_eal_max = extrema(relative_eal[1:end-1])
figure();
plot(r_surf,intDD_θ/intDD_θ[end]);
plot(r_surf,exp.(-(μ0.-r_surf)/λe));
plot(r_surf,exp.(-(μ0.-r_surf)/(r_eal_min*λe)),color="red")
plot(r_surf,exp.(-(μ0.-r_surf)/(r_eal_max*λe)),color="red")
legend(["cylinder","plane \$\\lambda_e\$","plane \$r_{\\max}\\lambda_e\$","plane \$r_{\\min}\\lambda_e\$"])
xlabel("raduis [\$\\mu\$m]")
ylabel("relative gain [\$\\mu\$m]")
savefig("plans_vs_cylinder.png")
savefig("plans_vs_cylinder.pdf")


fig = figure(figsize=[5,5])
ax = fig.add_axes([0.1,0.1,0.8,0.8],polar=true)
ax.pcolormesh(θ,r_surf,intDD,edgecolors="face")



extrema(DD)
extrema(DD_simple)
