figure(figsize=[12, 10]); 

for (i_plot,i_data) in zip(1:4,Int64.(round.(collect(LinRange(1,Ndata,4)))))
    local ax1 = subplot(2,2,i_plot)
    title(string("Eph =", Int64(round(hν[i_data]))," [eV]"),fontsize=14)
    symbol_h = Symbol(string("hν_",Int64(round(hν[i_data])))) # mySym[i] # :hν_365
    plot(dictAllData[symbol_h][1][plot_sym],dictAllData[symbol_h][1].Sbg,label="background"); 
    plot(dictAllData[symbol_h][1][plot_sym],dictAllData[symbol_h][1].Sbg+dictAllData[symbol_h][1].SpectrumA_1,label="noise free spectrum"); 
    # scatter(dictAllData[symbol_h][1][plot_sym],dictAllData[symbol_h][1].Snoisy,label="noisy spectrum"); 
    plot(dictAllData[symbol_h][1][plot_sym],dictAllData[symbol_h][1].Snoisy,label="noisy spectrum"); 
    xlabel("kinetic energy [eV]",fontsize=14); ylabel("spectrum [a.u.]",fontsize=14) 
    xticks(fontsize=14); yticks(fontsize=14); 
    legend(fontsize=14) # ,bbox_to_anchor=[0.6,0.95],loc=2,borderaxespad=0
    if (plot_sym==:Be)
        xlabel("binding energy [eV]",fontsize=14); 
        ax1.invert_xaxis();
    end
    ax1.ticklabel_format(axis="y",style="sci",scilimits=(-1,1),useOffset=true)
    ax1.yaxis.offsetText.set_size(14)
    ax1.xaxis.offsetText.set_size(14)
    ax1.tick_params(labelsize=14)
end

tight_layout(pad=1.0, w_pad=0.5, h_pad=0.2)

if SAVE_FIG
    savefig(string(save_folder,exp_tag,"/full_measurement_model.png"))
    savefig(string(save_folder,exp_tag,"/full_measurement_model.pdf"))
end