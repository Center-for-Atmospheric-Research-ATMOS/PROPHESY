# square root matrix of the generative covariance matrix (the covariance in the distribution used for generating new samples)
σw = 1.0e-3 # 0.001; # small compared with the amplitude of the state 
w = σw*ones(Cdouble,N_lowres); # not optimal because we know that the concentration varies more in the region near the surface rather than deep in the sample
Γsqrt = real(sqrt(corrCovariance(w;cor_len=10.0)));
p0 = 0.5 # 8.0*0.05 # 0.02; #starting acceptance rate of uphill moves
Ns = 1000000;

ΓIinv = zeros(Cdouble,Ndata,Ndata,Nnoise);
for i in 1:Nnoise
    [ΓIinv[k,k,i] = 1.0/(ΓI[k,k,i]+(σB*σεH[k])^2) for k in 1:Ndata]
end

μρ_I = zeros(Cdouble,N_lowres,Nnoise);
Γρ_I = zeros(Cdouble,N_lowres,N_lowres,Nnoise);
μρ_I_un = zeros(Cdouble,N_lowres,Nnoise);
Γρ_I_un = zeros(Cdouble,N_lowres,N_lowres,Nnoise);
for i in 1:Nnoise
    # conditional to data and model
    local ρ_all = samplePosterior(μρ[2:N0_lowres],Γsqrt,p0*ones(Cdouble,Ns),y_tilde[i,:],yd,ΓIinv[:,:,i],Γd_lowres,H_tilde,D_tilde;Ns=Ns);

    # error marginalization
    global ΓHΓyinv .= 0.0;
    for k in 1:Ndata
        # ΓHΓyinv = ΓHΓyinv + diagm(diag(ΓH[2:N0_lowres,2:N0_lowres,k])/ΓI[k,k,i])
        ΓHΓyinv = ΓHΓyinv + ΓH[2:N0_lowres,2:N0_lowres,k]/(ΓI[k,k,i]+(σB*σεH[k])^2)
    end
    local ρ_all_un = samplePosteriorMargin(μρ_un[2:N0_lowres],Γsqrt,p0*ones(Cdouble,Ns),y_tildeμ[i,:],yd,ΓIinv[:,:,i],Γd_lowres,μH_tilde,D_tilde,ΓHΓyinv;Ns=Ns);

    # compute a covariance matrix from the samples 
    μρ_I[:,i] = dropdims(mean(ρ_all,dims=1),dims=1);
    Γρ_I[:,:,i] = cov(ρ_all);
    μρ_I_un[:,i] = dropdims(mean(ρ_all_un,dims=1),dims=1);
    Γρ_I_un[:,:,i] = cov(ρ_all_un);
end

μμρ_I = dropdims(mean(μρ_I,dims=2),dims=2);
μΓρ_I = dropdims(mean(Γρ_I,dims=3),dims=3);

μμρ_I_un = dropdims(mean(μρ_I_un,dims=2),dims=2);
μΓρ_I_un = dropdims(mean(Γρ_I_un,dims=3),dims=3);

figure()
plot(1000.0(r_lowres.-μ0),reverse(μρ),color="tab:red",label="CP estimation")
fill_between(1000.0(r_lowres[end-N0_lowres+1:end-1].-μ0),reverse(μρ[2:N0_lowres]-sqrt.(diag(μΓρ_I))),reverse(μρ[2:N0_lowres]+sqrt.(diag(μΓρ_I))),alpha=0.5,color="tab:red",label="uncertainty posterior")
plot(1000.0(r_lowres[end-N0_lowres+1:end-1].-μ0),reverse(μμρ_I),color="tab:blue",label="sampling mean")
fill_between(1000.0(r_lowres[end-N0_lowres+1:end-1].-μ0),reverse(μμρ_I-sqrt.(diag(μΓρ_I))),reverse(μμρ_I+sqrt.(diag(μΓρ_I))),alpha=0.5,color="tab:blue",label="uncertainty posterior")
plot(1000.0(r.-μ0),reverse(ρA_1),color="tab:green",label="GT")
legend()
legend(fontsize=14)
xlabel("depth [nm]",fontsize=14)
xticks(fontsize=14)
ylabel("concentration [a.u.]",fontsize=14)
yticks(fontsize=14)




figure()
plot(1000.0(r_lowres.-μ0),reverse(μρ_un),color="tab:red",label="CP estimation")
fill_between(1000.0(r_lowres[end-N0_lowres+1:end-1].-μ0),reverse(μρ_un[2:N0_lowres]-sqrt.(diag(μΓρ_I_un))),reverse(μρ_un[2:N0_lowres]+sqrt.(diag(μΓρ_I_un))),alpha=0.5,color="tab:red",label="uncertainty posterior")

plot(1000.0(r_lowres[end-N0_lowres+1:end-1].-μ0),reverse(μμρ_I_un),color="tab:blue",label="sampling mean")
fill_between(1000.0(r_lowres[end-N0_lowres+1:end-1].-μ0),reverse(μμρ_I_un-sqrt.(diag(μΓρ_I_un))),reverse(μμρ_I_un+sqrt.(diag(μΓρ_I_un))),alpha=0.5,color="tab:blue",label="uncertainty posterior")

plot(1000.0(r.-μ0),reverse(ρA_1),color="tab:green",label="GT")
legend()
legend(fontsize=14)
xlabel("depth [nm]",fontsize=14)
xticks(fontsize=14)
ylabel("concentration [a.u.]",fontsize=14)
yticks(fontsize=14)


# fill_between(1000.0(r_lowres.-μ0),reverse(μρ_un-sqrt.(diag(Γρ_un))),reverse(μρ_un+sqrt.(diag(Γρ_un))),alpha=0.5,color="tab:orange",label="uncertainty marginal")
