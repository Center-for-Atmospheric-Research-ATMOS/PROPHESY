## compute distances
# DD_sphere_all      = d_sphere_P(r,θ,φ,x0_near,y0_near,z0_near,μ0);
# DD_sphere_away_all = d_sphere_P(r,θ,φ,x0_far,y0_far,z0_far,μ0);
d_offset = (1.1μ0,0.0*1.1μ0,1.1μ0);
DD_sphere_all      = d_sphere_P(r,θ,φ,d_offset[1],d_offset[2],d_offset[3],μ0);
DD_sphere_away_all = d_sphere_P(r,θ,φ,10.0μ0,0.0,0.0,μ0);
DD_sphere_90      = DD_sphere_all[:,:,25]; # DD_sphere_all[:,:,128];
DD_sphere_away_90 = DD_sphere_away_all[:,:,25]; # DD_sphere_away_all[:,:,128];
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



# plot a 3D surface with face color

# xx = cos.(θ)*(sin.(φ)');
# yy = sin.(θ)*(sin.(φ)');
# zz = ones(Cdouble,length(θ))*(cos.(φ)');
idx_r = 40;
maxDist = maximum(DD_sphere_all);
fig,ax,myCm = plot_sphere(r[idx_r],θ,φ,DD_sphere_all[idx_r,:,:]./maxDist;fig=963,cm_heat=PyPlot.cm.YlGnBu_r)

idx_r = 50
fig = figure(962,figsize=[9,5])
axx = fig.add_subplot(121, projection="3d");
axmp1,myCm = plot_sphere_ax(axx,r[idx_r],θ,φ,DD_sphere_all[idx_r,:,:]./maxDist);
axmp2,myCm = plot_sphere_ax(axx,1.0,θ,φ,0.0DD_sphere_all[idx_r,:,:]./maxDist;d_offset=d_offset);
cax = fig.add_axes([0.01,0.25,0.03,0.3])
cb = fig.colorbar(axmp1, orientation="vertical", cax=cax)
cb.update_ticks()

ayy = fig.add_subplot(122, projection="3d");
idx_φ = 14; # 9; #17; # 26;
axmp3,myCm = plot_cone_ax(ayy,r,θ,φ[idx_φ],DD_sphere_all[:,:,idx_φ]./maxDist)
idx_φ = 37; # 43; #17; # 26;
axmp4,myCm = plot_cone_ax(ayy,r,θ,φ[idx_φ],DD_sphere_all[:,:,idx_φ]./maxDist)
axmp5,myCm = plot_sphere_ax(ayy,1.0,θ,φ,0.0DD_sphere_all[idx_r,:,:]./maxDist;d_offset=d_offset);
xlim(-r[end],r[end])
ylim(-r[end],r[end])
zlim(-r[end],r[end])
cax2 = fig.add_axes([0.51,0.25,0.03,0.3])
cb2 = fig.colorbar(axmp3, orientation="vertical", cax=cax2)
cb2.update_ticks()


tight_layout(pad=1.0, w_pad=2.0, h_pad=0.2)
axx.text2D(0.05, 0.9, "a)", transform=axx.transAxes,fontsize=14)
ayy.text2D(1.05, 0.9, "b)", transform=axx.transAxes,fontsize=14)

# axx.annotate("a)", xy=(3, 1),  xycoords="data", xytext=(0.0, 1.0), textcoords="axes fraction", color="black",fontsize=14)
# ayy.annotate("b)", xy=(3, 1),  xycoords="data", xytext=(0.0, 1.0), textcoords="axes fraction", color="black",fontsize=14)

# savefig("distance_sphere.png")
# savefig("distance_sphere.pdf")
