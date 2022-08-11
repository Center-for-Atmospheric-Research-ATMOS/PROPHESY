figure(figsize=[12, 10]); 
for (i_plot,i_data) in zip(1:4,Int64.(round.(collect(LinRange(1,Ndata,4)))))
    local ax1 = subplot(2,2,i_plot)
    local symbol_h = key_symbol[i_data]
    title(string("Eph =", Int64(round(dictAllData[symbol_h][1][:hν][1]))," [eV]"),fontsize=14)
    plot(dictAllData[symbol_h][1][plot_sym],Sbg_est[symbol_h],color=color_array[1],label="estimated background"); 
    plot(dictAllData[symbol_h][1][plot_sym],dictAllData[symbol_h][1].Sbg,color=color_array[2],label="background gt"); 
    plot(dictAllData[symbol_h][1][plot_sym],dictAllData[symbol_h][1].SpectrumA_1,color=color_array[3],label="SOI gt"); 
    scatter(dictAllData[symbol_h][1][plot_sym],dictAllData[symbol_h][1].Snoisy-Sbg_est[symbol_h],color=color_array[1],label="estimated SOI"); 
    xlabel("kinetic energy [eV]",fontsize=14); ylabel("spectrum [count]",fontsize=14) 
    xticks(fontsize=14); yticks(fontsize=14); 
    legend(fontsize=14)
    if (plot_sym==:Be)
        xlabel("binding energy [eV]",fontsize=14); 
        #ax = gca()
        ax1.invert_xaxis();
    end
end
tight_layout(pad=1.0, w_pad=0.5, h_pad=0.2)
# savefig(string(data_tag,"_background_removal.pdf"))
# savefig(string(data_tag,"_background_removal.png"))

