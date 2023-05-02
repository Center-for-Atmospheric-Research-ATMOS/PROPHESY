figure(figsize=[12, 10]); 
# for j in  1:Ndata # 1:5:Ndata
for (i_plot,i_data) in zip(1:4,Int64.(round.(collect(LinRange(1,Ndata,4)))))
    local μKe0=50.0
    local μKe1=1200.0
    local symbol_h = key_symbol[i_data]; 
    local be = dictAllData[symbol_h][1][plot_sym] #dictAllData[symbol_h][1][:Be];
    
    local μKe = dictAllData[symbol_h][1][:μKe];
    # partial cross section (one for each chemical state)
    local σ_peak_1 = (1.0/sqrt(2.0π*σ_be[1]^2))*exp.(-(be.-μBe[1]).^2/(2.0σ_be[1]^2));
    local σ_peak_2 = (1.0/sqrt(2.0π*σ_be[2]^2))*exp.(-(be.-μBe[2]).^2/(2.0σ_be[2]^2));
    local σ_peak_3 = (1.0/sqrt(2.0π*σ_be[3]^2))*exp.(-(be.-μBe[3]).^2/(2.0σ_be[3]^2));
    # quantity of chemical states
    p1 = 0.85 .+ (0.77-0.85)*(μKe[1].-μKe0)./(μKe1-μKe0);
    p2 = 0.125 .+ (0.12-0.125)*(μKe[1].-μKe0)./(μKe1-μKe0);
    p3 = 1.0-(p1+p2);

    # plotting
    ax1 = subplot(2,2,i_plot)
    title(string("Eph =", Int64(round(dictAllData[symbol_h][1][:hν][1]))," [eV]"),fontsize=14)
    plot(be,p1*σ_peak_1+p2*σ_peak_2+p3*σ_peak_3,color=color_array[3],label="GT") # i_plot
    plot(be,S_cs_dens[symbol_h],color=color_array[1],label="estimation") # i_plot
    scatter(be,(dictAllData[symbol_h][1][:Snoisy]-dictAllData[symbol_h][1][:Sbg])/(dKe*sum(dictAllData[symbol_h][1][:Snoisy]-dictAllData[symbol_h][1][:Sbg])),color=color_array[2],label="scaled data") # i_plot

    xlabel("kinetic energy [eV]",fontsize=14); 
    ylabel("cross section density [eV\$^{-1}\$]",fontsize=14) 
    xticks(fontsize=14); yticks(fontsize=14); 
    legend(fontsize=14,bbox_to_anchor=[0.05,0.95],loc=2,borderaxespad=0) #

    if (plot_sym==:Be)
        xlabel("binding energy [eV]",fontsize=14); 
        ax1.invert_xaxis();
    end
end
tight_layout(pad=1.0, w_pad=0.5, h_pad=0.2)
# savefig(string(data_tag,"_cross_section_density.pdf"))
# savefig(string(data_tag,"_cross_section_density.png"))

# TODO: plot the variation of the peak estimates w.r.t. the photon energy (maybe use the case with 20 measurements to make the trend clear... it's getting worse with the probing depth)

figure()
scatter(hν_list,μt[:,1].-μBe[1],label="E\$_b\$=290.2 [eV]")
scatter(hν_list,μt[:,2].-μBe[2],label="E\$_b\$=292.0 [eV]")
scatter(hν_list,μt[:,3].-μBe[3],label="E\$_b\$=293.0 [eV]")
ylim(-0.1,0.1)
xlim(250.0,2000.0)
ax = gca()
ax.ticklabel_format(axis="y",style="sci",scilimits=(-1,2),useOffset=true)
ax.yaxis.offsetText.set_size(14)
ax.xaxis.offsetText.set_size(14)
legend(fontsize=14,bbox_to_anchor=[0.6,0.95],loc=2,borderaxespad=0)
xticks(fontsize=14); yticks(fontsize=14); 
ax.tick_params(labelsize=14)
xlabel("photon energy [eV]",fontsize=14)
ylabel("binding energy deviation [eV]",fontsize=14)
tight_layout(pad=1.0, w_pad=0.5, h_pad=0.2)
# savefig(string(data_tag,"_binding_energy.pdf"))
# savefig(string(data_tag,"_binding_energy.png"))
