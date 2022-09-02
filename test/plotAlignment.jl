# compute the fit for the O1s data
global sum_xixi = 0.0; 
global sum_xi = 0.0; 
global sum_yixi = 0.0; 
global sum_yi = 0.0; 
global sum_1 = 0.0;
for (key_O1s,_) in α_noiseO1s
    local xi = log.(α_noiseO1s[key_O1s])
    local yi = log.(α_ratioO1s[key_O1s])
    global sum_xixi = sum_xixi + sum(xi.^2)
    global sum_xi   = sum_xi + sum(xi)
    global sum_yi   = sum_yi + sum(yi)
    global sum_yixi = sum_yixi + sum(xi.*yi)
    global sum_1    = sum_1  + length(xi)
end
for (key_O1s,_) in α_noiseO1sS2p
    local xi = log.(α_noiseO1sS2p[key_O1s])
    local yi = log.(α_ratioO1sS2p[key_O1s])
    global sum_xixi = sum_xixi + sum(xi.^2)
    global sum_xi   = sum_xi + sum(xi)
    global sum_yi   = sum_yi + sum(yi)
    global sum_yixi = sum_yixi + sum(xi.*yi)
    global sum_1    = sum_1  + length(xi)
end
μab =  (inv([sum_xixi sum_xi; sum_xi sum_1])*[sum_yixi;sum_yi])';
α_noise_log_fit = log.(10.0*[1.0e-10 2.0e-10 5.0e-10 1.0e-9 2.0e-9 5.0e-9 1.0e-8])
α_ratio_log_fit = dropdims(μab*[α_noise_log_fit; ones(Cdouble,length(α_noise_log_fit))'],dims=1)

figure(figsize=[12, 10])
ax1 = subplot(221)
for i in 1:length(α_noiseO1s)
    local plot_sym = Symbol(data_filesO1s[i][1:end-5]);
    scatter(1.0e8mean(α_noiseO1s[plot_sym]),mean(α_ratioO1s[plot_sym]),s=100,color=color_array[i],alpha=0.5)
    scatter(1.0e8α_noiseO1s[plot_sym],α_ratioO1s[plot_sym],label=replace(data_filesO1s[i][1:end-5],"_"=>" ")) #.^2.5
end
plot(1.0e8exp.(dropdims(α_noise_log_fit,dims=1)),exp.(α_ratio_log_fit),label="O1s fit")
# xlim(0.1,10.0)
xlim(0.05,50.0)
ylim(0.1,10.0)
xlabel("model estimation [cm\$^{-2}\$]",fontsize=14); 
ylabel("liq O1s/gas O1s",fontsize=14) 
xticks(fontsize=14); yticks(fontsize=14); 
ax1.ticklabel_format(axis="y",style="sci",scilimits=(-1,1),useOffset=true)
ax1.yaxis.offsetText.set_size(14)
ax1.xaxis.offsetText.set_size(14)
xscale("log")
yscale("log")
# legend(fontsize=12)
legend(fontsize=12,borderpad=0.4,borderaxespad=0.2,handletextpad=0.5,handlelength=1.0,framealpha=0.4)

ax3 = subplot(223)
for i in 1:length(α_noiseO1sS2p)
    local plot_sym = Symbol(data_filesO1sS2p[i][1:end-5]);
    scatter(1.0e8mean(α_noiseO1sS2p[plot_sym]),mean(α_ratioO1sS2p[plot_sym]),s=100,color=color_array[i],alpha=0.5)
    scatter(1.0e8α_noiseO1sS2p[plot_sym],α_ratioO1sS2p[plot_sym],label=replace(data_filesO1sS2p[i][1:end-5],"_"=>" ")) #.^2.5
end
plot(1.0e8exp.(dropdims(α_noise_log_fit,dims=1)),exp.(α_ratio_log_fit),label="O1s fit")
xlim(0.05,50.0)
ylim(0.1,10.0)
xlabel("model estimation [cm\$^{-2}\$]",fontsize=14); 
ylabel("liq O1s/gas O1s",fontsize=14) 
xticks(fontsize=14); yticks(fontsize=14); 
ax3.ticklabel_format(axis="y",style="sci",scilimits=(-1,1),useOffset=true)
ax3.yaxis.offsetText.set_size(14)
ax3.xaxis.offsetText.set_size(14)
xscale("log")
yscale("log")
# legend(fontsize=12)
legend(fontsize=12,borderpad=0.4,borderaxespad=0.2,handletextpad=0.5,handlelength=1.0,framealpha=0.4)

ax2 = subplot(222)
for i in 1:length(α_noiseC1s)
    local plot_sym = Symbol(data_filesC1s[i][1:end-5]);
    scatter(1.0e8α_noiseC1s[plot_sym],α_ratioC1s[plot_sym],label=replace(data_filesC1s[i][1:end-5],"_"=>" "),color=color_array[i])
    scatter(1.0e8mean(α_noiseC1s[plot_sym]),mean(α_ratioC1s[plot_sym]),s=100,color=color_array[i],alpha=0.5) # .^2.5
end
plot(1.0e8exp.(dropdims(α_noise_log_fit,dims=1).-1.25),exp.(α_ratio_log_fit),label="O1s fit (scaled)")
xlim(0.01,10.0)
ylim(0.1,10.0)
xlabel("model estimation [cm\$^{-2}\$]",fontsize=14); 
ylabel("liq O1s/gas O1s",fontsize=14) 
xticks(fontsize=14); yticks(fontsize=14); 
ax2.ticklabel_format(axis="y",style="sci",scilimits=(-1,1),useOffset=true)
ax2.yaxis.offsetText.set_size(14)
ax2.xaxis.offsetText.set_size(14)
xscale("log")
yscale("log")
# legend(fontsize=12)
legend(fontsize=12,borderpad=0.4,borderaxespad=0.2,handletextpad=0.5,handlelength=1.0,framealpha=0.4)

ax4 = subplot(224)
for i in 1:length(α_noiseS2p)
    local plot_sym = Symbol(data_filesS2p[i][1:end-5]);
    scatter(1.0e8α_noiseS2p[plot_sym],α_ratioS2p[plot_sym],label=replace(data_filesS2p[i][1:end-5],"_"=>" "),color=color_array[i])
    scatter(1.0e8mean(α_noiseS2p[plot_sym]),mean(α_ratioS2p[plot_sym]),s=100,color=color_array[i],alpha=0.5) # .^2.5
end
plot(1.0e8exp.(dropdims(α_noise_log_fit.-5.75,dims=1)),exp.(α_ratio_log_fit),label="O1s fit (scaled)")
xlim(1.0e-4,0.1)
ylim(0.1,10.0)
xlabel("model estimation [cm\$^{-2}\$]",fontsize=14); 
ylabel("liq O1s/gas O1s",fontsize=14) 
xticks(fontsize=14); yticks(fontsize=14); 
ax4.ticklabel_format(axis="y",style="sci",scilimits=(-1,1),useOffset=true)
ax4.yaxis.offsetText.set_size(14)
ax4.xaxis.offsetText.set_size(14)
xscale("log")
yscale("log")
# legend(fontsize=12)
legend(fontsize=12,borderpad=0.4,borderaxespad=0.2,handletextpad=0.5,handlelength=1.0,framealpha=0.4)
tight_layout(pad=1.0, w_pad=0.2, h_pad=0.2)

ax1.text(-0.12, 0.9, "a)", transform=ax1.transAxes,fontsize=16)
ax2.text(-0.07, 0.9, "b)", transform=ax2.transAxes,fontsize=16)
ax3.text(-0.12, 0.9, "c)", transform=ax3.transAxes,fontsize=16)
ax4.text(-0.07, 0.9, "d)", transform=ax4.transAxes,fontsize=16)

if FLAG_SAVE_PLOT
    # savefig(string("../data/TK/","ratio_vs_model_O1s_C1s_S2p_mean.png"))
    # savefig(string("../data/TK/","ratio_vs_model_O1s_C1s_S2p_mean.pdf"))
    savefig(string("../data/TK/","ratio_vs_model_O1s_C1s_S2p_more_O1s_with_mean_units_new_cross_section.png"))
    savefig(string("../data/TK/","ratio_vs_model_O1s_C1s_S2p_more_O1s_with_mean_units_new_cross_section.pdf"))
end



    