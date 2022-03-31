## compute distances
DD = d_cylinder_P(r,θ,y,x0_near,y0_near,z0_near,μ0);
DD_simple = d_cylinder_P_simple(r,θ,y,x0_near,y0_near,z0_near,μ0);

DD_away = d_cylinder_P(r,θ,y,x0_far,y0_far,z0_far,μ0);
DD_simple_away = d_cylinder_P_simple(r,θ,y,x0_far,y0_far,z0_far,μ0);


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


## plot
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
