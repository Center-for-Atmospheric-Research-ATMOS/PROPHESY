J_max,J_min = extrema(J_all)
fig,ax,pcm = imshowData(10,1.0e3(μ0.-r),1.0e3(μ0.-r),J_all;_norm=:NoNorm,_vmin=J_min,_vmax=J_max,_edgecolors="face",_shading="None")
ax.ticklabel_format(axis="both",style="sci",scilimits=(-1,2),useOffset=true) #,useMathText=true
xlabel("surface distance [nm]",fontsize=14)
ylabel("surface distance [nm]",fontsize=14)
xticks(fontsize=14); 
yticks(fontsize=14);

rc("ytick",color="white")
cax   = fig.add_axes([0.8, .5, 0.03, 0.4])
cbar  = fig.colorbar(pcm, orientation="vertical", cax=cax, shrink=0.6)
cbar.formatter.set_powerlimits((-1,2))
cbar.update_ticks()
cbar.set_label("sensitivity [(m\$^{-3}\$)\$^{-2}\$]",fontsize=14, color="white")
cbar.ax.tick_params(labelsize=14)
cbar.formatter.set_powerlimits((-1,2))
cbar.ax.yaxis.offsetText.set_size(14)
cbar.ax.yaxis.set_tick_params(color="white")
cbar.outline.set_edgecolor("white")
rc("ytick",color="black")
tight_layout(pad=1.0, w_pad=0.5, h_pad=0.2)
# savefig(string(data_tag,"_sensitivity_matrix.pdf"))
# savefig(string(data_tag,"_sensitivity_matrix.png"))



# F_j_all  = svd(J_all);
# figure()
# plot(r,F_j_all.U[:,1])


# S_approx = zeros(Cdouble,Nr,Nr);
# S_approx[1,1] = F_j_all.S[1];
# J_approx = F_j_all.U*S_approx*F_j_all.Vt

# figure();
# imshow(J_all)
# colorbar()

# figure();
# imshow(J_approx)
# colorbar()


# figure();
# imshow(abs.(J_approx-J_all))
# colorbar()

# figure();
# imshow(abs.(F_j_all.U[:,1:Ndata-1]*F_j_all.Vt[1:Ndata-1,:]))
# colorbar()