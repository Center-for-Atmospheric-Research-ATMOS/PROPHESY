## compute distances
# DD_sphere_all      = d_sphere_P(r,θ,φ,x0_near,y0_near,z0_near,μ0);
# DD_sphere_away_all = d_sphere_P(r,θ,φ,x0_far,y0_far,z0_far,μ0);
DD_sphere_all      = d_sphere_P(r,θ,φ,1.1μ0,0.0,0.0,μ0);
DD_sphere_away_all = d_sphere_P(r,θ,φ,10.0μ0,0.0,0.0,μ0);
DD_sphere_90      = DD_sphere_all[:,:,128];
DD_sphere_away_90 = DD_sphere_away_all[:,:,128];
DD_sphere_0      = DD_sphere_all[:,:,1];
DD_sphere_away_0 = DD_sphere_away_all[:,:,1];

## plot
fig = figure(figsize=[9,5])
ax1 = subplot(121,polar=true)
ax1.set_rticks([μ0/4, 2μ0/4, 3μ0/4, μ0])
yticks(fontsize=12)
xticks(fontsize=12)
ax1.set_rlabel_position(-22.5)
pcm1 = ax1.pcolormesh(θ,r,DD_sphere_0,edgecolors="face")
cax1 = fig.add_axes([0.08, .35, 0.02, 0.3])
cb1 = fig.colorbar(pcm1, orientation="vertical", cax=cax1, shrink=0.6)
cb1.set_label("distance [\$\\mu\$m]", fontsize=12) # , color="white"

ax2 = subplot(122,polar=true)
ax2.set_rticks([μ0/4, 2μ0/4, 3μ0/4, μ0])
yticks(fontsize=12)
xticks(fontsize=12)
ax2.set_rlabel_position(-22.5)
pcm2 = ax2.pcolormesh(θ,r,DD_sphere_90,edgecolors="face")
cax2 = fig.add_axes([0.58, .35, 0.02, 0.3])
cb2 = fig.colorbar(pcm2, orientation="vertical", cax=cax2, shrink=0.6)
cb2.set_label("distance [\$\\mu\$m]", fontsize=12)


tight_layout(pad=1.0, w_pad=2.0, h_pad=0.2)

ax1.annotate("a)", xy=(3, 1),  xycoords="data", xytext=(0.0, 1.0), textcoords="axes fraction", color="black",fontsize=14)
ax2.annotate("b)", xy=(3, 1),  xycoords="data", xytext=(0.0, 1.0), textcoords="axes fraction", color="black",fontsize=14)

# savefig("distance_cylinder_near.png")
# savefig("distance_cylinder_near.pdf")


## plot
fig = figure(figsize=[9,5])
ax1 = subplot(121,polar=true)
pcm1 = ax1.pcolormesh(θ,r,DD_sphere_away_0,edgecolors="face")
cax1 = fig.add_axes([0.08, .35, 0.02, 0.3])
cb1 = fig.colorbar(pcm1, orientation="vertical", cax=cax1, shrink=0.6)
cb1.set_label("distance [\$\\mu\$m]", fontsize=10)


ax2 = subplot(122,polar=true)
pcm2 = ax2.pcolormesh(θ,r,DD_sphere_away_90,edgecolors="face")
cax2 = fig.add_axes([0.58, .35, 0.02, 0.3])
cb2 = fig.colorbar(pcm2, orientation="vertical", cax=cax2, shrink=0.6)
cb2.set_label("distance [\$\\mu\$m]", fontsize=10)
tight_layout(pad=1.0, w_pad=2.0, h_pad=0.2)

ax1.annotate("a)", xy=(3, 1),  xycoords="data", xytext=(0.0, 1.0), textcoords="axes fraction", color="black",fontsize=14)
ax2.annotate("b)", xy=(3, 1),  xycoords="data", xytext=(0.0, 1.0), textcoords="axes fraction", color="black",fontsize=14)

# savefig("distance_cylinder_far.png")
# savefig("distance_cylinder_far.pdf")



# n = 300;
# θ = LinRange(0.0,π,n);
# ϕ = LinRange(0.0,2π,n);
# xx = @. sin(θ)*cos(φ');
# yy = @. sin(θ)*sin(φ');
# zz = @. cos(θ)*ones(K)';

xx = @. sin(φ)*cos(θ');
yy = @. sin(φ)*sin(θ');
zz = @. cos(φ)*ones(J)';

maxDist = maximum(DD_sphere_all);
myFaceColor = PyPlot.cm.YlGnBu_r(DD_sphere_all[25,:,:]./maxDist);
figure(); plot_surface(xx,yy,zz,cmap=PyPlot.cm.YlGnBu_r, facecolors=myFaceColor)

