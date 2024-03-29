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
α_noise_log_fit = log.(200.0*[1.0e-10 2.0e-10 5.0e-10 1.0e-9 2.0e-9 5.0e-9 1.0e-8])
α_ratio_log_fit = dropdims(μab*[α_noise_log_fit; ones(Cdouble,length(α_noise_log_fit))'],dims=1)

TWO_COLUMN = true
if TWO_COLUMN
    FONTSIZE = 16
else
    FONTSIZE = 14
end
if TWO_COLUMN
    LEGEND_FONTSIZE = 12
else
    LEGEND_FONTSIZE = 12
end
LINEWIDTH = 2.5
# figure(figsize=[12, 10])
if TWO_COLUMN
    figure(figsize=[12, 8])
else
    figure(figsize=[12, 10])
end
ax1 = subplot(221)
# for i in 1:length(α_noiseO1s)
idx_color = 1
for i in [1; 3; 2]
    local plot_sym = Symbol(data_filesO1s[i][1:end-5]);
    scatter(1.0e8mean(α_noiseO1s[plot_sym]),mean(α_ratioO1s[plot_sym]),s=100,color=color_array[idx_color],alpha=0.5)
    scatter(1.0e8α_noiseO1s[plot_sym],α_ratioO1s[plot_sym],label=replace(data_filesO1s[i][1:end-5],"_"=>" ")[5:end],color=color_array[idx_color]) #.^2.5
    global idx_color = idx_color + 1
end
α_noise_log_fit = log.([5e-8 0.9e-5])
α_ratio_log_fit = dropdims(μab*[α_noise_log_fit; ones(Cdouble,length(α_noise_log_fit))'],dims=1)
plot(1.0e8exp.(dropdims(α_noise_log_fit,dims=1)),exp.(α_ratio_log_fit),label="O1s fit")
# xlim(0.1,10.0)
# xlim(0.05,50.0)
xlim(1.0,1000.0)
ylim(0.1,10.0)
# xlabel("model estimation [cm\$^{-2}\$ s]",fontsize=FONTSIZE); 
# ylabel("liq O1s/gas O1s",fontsize=FONTSIZE) 
ylabel("LGPAR",fontsize=FONTSIZE) 
xticks(fontsize=FONTSIZE); yticks(fontsize=FONTSIZE); 
ax1.ticklabel_format(axis="y",style="sci",scilimits=(-1,1),useOffset=true)
ax1.yaxis.offsetText.set_size(FONTSIZE)
ax1.xaxis.offsetText.set_size(FONTSIZE)
xscale("log")
yscale("log")
# legend(fontsize=12)
legend(fontsize=LEGEND_FONTSIZE,borderpad=0.4,borderaxespad=0.2,handletextpad=0.5,handlelength=1.0,framealpha=0.3)

ax3 = subplot(223)
# for i in 1:length(α_noiseO1sS2p)
idx_color = 1
for i in [5; 4; 2; 3; 1]
    local plot_sym = Symbol(data_filesO1sS2p[i][1:end-5]);
    scatter(1.0e8mean(α_noiseO1sS2p[plot_sym]),mean(α_ratioO1sS2p[plot_sym]),s=100,color=color_array[idx_color],alpha=0.5)
    scatter(1.0e8α_noiseO1sS2p[plot_sym],α_ratioO1sS2p[plot_sym],label=replace(data_filesO1sS2p[i][1:end-5],"_"=>" ")[5:end-8],color=color_array[idx_color]) #.^2.5
    global idx_color = idx_color + 1
end
α_noise_log_fit = log.([2e-8 2.0e-6])
α_ratio_log_fit = dropdims(μab*[α_noise_log_fit; ones(Cdouble,length(α_noise_log_fit))'],dims=1)
plot(1.0e8exp.(dropdims(α_noise_log_fit,dims=1)),exp.(α_ratio_log_fit),label="O1s fit")
# xlim(0.05,50.0)
xlim(1.0,1000.0)
ylim(0.1,10.0)
xlabel("model estimation [cm\$^{-2}\$ s]",fontsize=FONTSIZE); 
# ylabel("liq O1s/gas O1s",fontsize=FONTSIZE) 
ylabel("LGPAR",fontsize=FONTSIZE) 
xticks(fontsize=FONTSIZE); yticks(fontsize=FONTSIZE); 
ax3.ticklabel_format(axis="y",style="sci",scilimits=(-1,1),useOffset=true)
ax3.yaxis.offsetText.set_size(14)
ax3.xaxis.offsetText.set_size(14)
xscale("log")
yscale("log")
# legend(fontsize=12)
legend(fontsize=LEGEND_FONTSIZE,borderpad=0.4,borderaxespad=0.2,handletextpad=0.5,handlelength=1.0,framealpha=0.3)

ax2 = subplot(222)
# for i in 1:length(α_noiseC1s)
idx_color = 2
for i in [4; 2; 3; 1]
    local plot_sym = Symbol(data_filesC1s[i][1:end-5]);
    scatter(1.0e8α_noiseC1s[plot_sym],α_ratioC1s[plot_sym],label=replace(data_filesC1s[i][1:end-5],"_"=>" ")[5:end],color=color_array[idx_color])
    scatter(1.0e8mean(α_noiseC1s[plot_sym]),mean(α_ratioC1s[plot_sym]),s=100,color=color_array[idx_color],alpha=0.5) # .^2.5
    global idx_color = idx_color + 1
end
α_noise_log_fit = log.([5e-8 1.0e-5])
α_ratio_log_fit = dropdims(μab*[α_noise_log_fit; ones(Cdouble,length(α_noise_log_fit))'],dims=1)
plot(1.0e8exp.(dropdims(α_noise_log_fit,dims=1).-1.25),exp.(α_ratio_log_fit),label="O1s fit (scaled)")
# xlim(0.01,10.0)
xlim(1.0,1000.0)
ylim(0.1,10.0)
# xlabel("model estimation [cm\$^{-2}\$ s]",fontsize=FONTSIZE); 
# ylabel("liq O1s/gas O1s",fontsize=FONTSIZE) 
# ylabel("LGPAR",fontsize=FONTSIZE) 
xticks(fontsize=FONTSIZE); yticks(fontsize=FONTSIZE); 
ax2.ticklabel_format(axis="y",style="sci",scilimits=(-1,1),useOffset=true)
ax2.yaxis.offsetText.set_size(14)
ax2.xaxis.offsetText.set_size(14)
xscale("log")
yscale("log")
# legend(fontsize=12)
legend(fontsize=LEGEND_FONTSIZE,borderpad=0.4,borderaxespad=0.2,handletextpad=0.5,handlelength=1.0,framealpha=0.3)

ax4 = subplot(224)
# for i in 1:length(α_noiseS2p)
idx_color = 1
for i in [2; 5; 3; 4; 1]
    local plot_sym = Symbol(data_filesS2p[i][1:end-5]);
    # scatter(1.0e8α_noiseS2p[plot_sym],α_ratioS2p[plot_sym],label=replace(data_filesS2p[i][1:end-5],"_"=>" "),color=color_array[i])
    scatter(1.0e8α_noiseS2p[plot_sym],α_ratioS2p[plot_sym],label=replace(data_filesS2p[i][1:end-5],"_"=>" ")[5:end],color=color_array[idx_color])
    scatter(1.0e8mean(α_noiseS2p[plot_sym]),mean(α_ratioS2p[plot_sym]),s=100,color=color_array[idx_color],alpha=0.5) # .^2.5
    global idx_color = idx_color + 1
end
α_noise_log_fit = log.([5e-8 2.0e-6])
α_ratio_log_fit = dropdims(μab*[α_noise_log_fit; ones(Cdouble,length(α_noise_log_fit))'],dims=1)
plot(1.0e8exp.(dropdims(α_noise_log_fit.-5.75,dims=1)),exp.(α_ratio_log_fit),label="O1s fit (scaled)")
# xlim(1.0e-4,0.1)
xlim(5.0e-4,1.0)
ylim(0.1,10.0)
xlabel("model estimation [cm\$^{-2}\$ s]",fontsize=FONTSIZE); 
# ylabel("liq O1s/gas O1s",fontsize=FONTSIZE) 
# ylabel("LGPAR",fontsize=FONTSIZE) 
xticks(fontsize=FONTSIZE); yticks(fontsize=FONTSIZE); 
ax4.ticklabel_format(axis="y",style="sci",scilimits=(-1,1),useOffset=true)
ax4.yaxis.offsetText.set_size(FONTSIZE)
ax4.xaxis.offsetText.set_size(FONTSIZE)
xscale("log")
yscale("log")
# legend(fontsize=12)
legend(fontsize=LEGEND_FONTSIZE,borderpad=0.4,borderaxespad=0.2,handletextpad=0.5,handlelength=1.0,framealpha=0.3)
tight_layout(pad=0.5, w_pad=0.2, h_pad=0.2)

ax1.text(-0.12, 0.9, "a)", transform=ax1.transAxes,fontsize=16)
ax2.text(-0.07, 0.9, "b)", transform=ax2.transAxes,fontsize=16)
ax3.text(-0.12, 0.9, "c)", transform=ax3.transAxes,fontsize=16)
ax4.text(-0.07, 0.9, "d)", transform=ax4.transAxes,fontsize=16)

ax1.text(0.65, 0.25, "O 1s", transform=ax1.transAxes,fontsize=16)
# ax2.text(0.65, 0.25, "C 1s", transform=ax2.transAxes,fontsize=16)
ax2.text(0.65, 0.5, "C 1s", transform=ax2.transAxes,fontsize=16)
ax3.text(0.65, 0.25, "O 1s for S 2p", transform=ax3.transAxes,fontsize=16)
ax4.text(0.65, 0.25, "S 2p", transform=ax4.transAxes,fontsize=16)

if FLAG_SAVE_PLOT
    # savefig(string("../data/TK/","ratio_vs_model_O1s_C1s_S2p_mean.png"))
    # savefig(string("../data/TK/","ratio_vs_model_O1s_C1s_S2p_mean.pdf"))
    # savefig(string("../data/TK/","ratio_vs_model_O1s_C1s_S2p_more_O1s_with_mean_units_new_cross_section.png"))
    # savefig(string("../data/TK/","ratio_vs_model_O1s_C1s_S2p_more_O1s_with_mean_units_new_cross_section.pdf"))
    # savefig(string("../data/TK/","two_col_ratio_vs_model_O1s_C1s_S2p_more_O1s_with_mean_units_new_cross_section.png"))
    # savefig(string("../data/TK/","two_col_ratio_vs_model_O1s_C1s_S2p_more_O1s_with_mean_units_new_cross_section.pdf"))
    # savefig(string("../data/TK/","two_col_ratio_vs_model_O1s_C1s_S2p_more_O1s_with_mean_units_new_cross_section_new_color_because_it_was_supposedly_confusing.png"))
    # savefig(string("../data/TK/","two_col_ratio_vs_model_O1s_C1s_S2p_more_O1s_with_mean_units_new_cross_section_new_color_because_it_was_supposedly_confusing.pdf"))
    savefig(string("../../../data/TK/","two_col_ratio_vs_model_O1s_C1s_S2p_more_O1s_with_mean_units_new_cross_section_new_color_because_it_was_supposedly_confusing_proper_unit.png"))
    savefig(string("../../../data/TK/","two_col_ratio_vs_model_O1s_C1s_S2p_more_O1s_with_mean_units_new_cross_section_new_color_because_it_was_supposedly_confusing_proper_unit.pdf"))
end



    