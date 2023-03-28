## load the packages used in the estimation
# plotting
using PyPlot
rc("text", usetex=true)
rc("figure",max_open_warning=50)
using myPlot

# data manipulation (loading, writing, etc)
using Printf

# scientific package from the official Julia repositories
using LinearAlgebra
using StatsBase

# modeling XPS
using XPSpack
using XPSinv
using XPSsampling


# create some data with a fairly fine discretization and another operator with a different resolution
include("dataGenCylinder.jl")


##TODO: compute the correction with the other measurement model (different discretization and slightly different penetration depth (maybe add a random element in lambda?))
# data correction due to the signal coming from the deep strates of the surface for which we will assume that the bulk value has been reach for the concentration
Δy = H_better[:,N0+1:end]*ρA_1[N0+1:end];
# NOTE: maybe this part should be discussed (what should be the value at the "surface"?)
δy = H_better[:,1]*ρA_1[1];

figure(); plot(y_data'); plot((y_data.-Δy')'); plot((y_data.-(Δy+δy)')')
plot(Δy); plot(δy)

y_tilde = y_data.-(Δy+δy)';
H_tilde = H_better[:,2:N0]


N=N0-1;
N_lowres = N0_lowres - 1

D_2nd = D2nd(N);
D_2nd_lowres = D2nd(N_lowres);


##
## other truncation
##


H0 = H_lowres[:,1];
μH0 = μH'[:,1];
H_tilde = H_lowres[:,2:N0_lowres];
μH_tilde = μH'[:,2:N0_lowres];
Hb = H_lowres[:,N0_lowres+1:end];
μHb = μH'[:,N0_lowres+1:end];
Δy = dropdims(sum(Hb,dims=2),dims=2)*ρA_1[end] # Hb*ρA_1[N0+1:end];
δy = H0*ρA_1[1];
Δyμ = dropdims(sum(μHb,dims=2),dims=2)*ρA_1[end] # Hb*ρA_1[N0+1:end];
δyμ = μH0*ρA_1[1];

y_tilde  = y_data.-(Δy+δy)';
y_tildeμ = y_data.-(Δyμ+δyμ)';

figure(); plot(y_data'); plot(y_tilde'); plot(y_tildeμ')

DN = D2nd(N_lowres+3); #+3
D0 = DN[:,1];
D_tilde = DN[:,2:N_lowres+1];
Db = DN[:,N_lowres+2:N_lowres+3];

Δyd = -dropdims(sum(Db,dims=2),dims=2)*ρA_1[end] # Db*ρA_1[N+2:N+3]
δyd = -D0*ρA_1[1];

yd = Δyd+δyd;

Htrunc = [H_tilde; D_tilde];

Htrunc_un = [μH_tilde; D_tilde; diagm(ones(Cdouble,N_lowres))]

ΓH0b = Array{Cdouble,3}(undef,Nr_lowres-N_lowres,Nr_lowres-N_lowres,Ndata);
for i in 1:Ndata
    ΓH0b[1,1,i] = ΓH[1,1,i]
    ΓH0b[1,2:end,i] = ΓH[1,N0_lowres+1:end,i]
    ΓH0b[2:end,1,i] = ΓH[N0_lowres+1:end,1,i]
    ΓH0b[2:end,2:end,i] = ΓH[N0_lowres+1:end,N0_lowres+1:end,i]
end
σεH = sqrt.(dropdims(sum(ΓH0b,dims=(1,2)),dims=(1,2)));
σB = 0.01 # 0.05;


# covariance matrix for the a priori distribution

Γprior = zeros(Cdouble,Nr,Nr);
Γprior_lowres = zeros(Cdouble,Nr_lowres,Nr_lowres);
cor_len = 5.0;
cor_len_lowres = cor_len/(Nr/Nr_lowres);
for i in 1:Nr
    # Γprior[i,i] = (1.0-ρA_1[i]+0.1)^2;
    Γprior[i,i] =  1.0
    for j in i+1:Nr
        Γprior[i,j] = Γprior[i,i]*exp(-(i-j)^2/(0.5*cor_len^2))
        Γprior[j,i] = Γprior[i,j]
    end
end
for i in 1:Nr_lowres
    # Γprior[i,i] = (1.0-ρA_1[i]+0.1)^2;
    Γprior_lowres[i,i] =  1.0
    for j in i+1:Nr_lowres
        Γprior_lowres[i,j] = Γprior_lowres[i,i]*exp(-(i-j)^2/(0.5*cor_len_lowres^2))
        Γprior_lowres[j,i] = Γprior_lowres[i,j]
    end
end
Γd = (0.01)^2*Γprior[2:N+2,2:N+2]; 
Γd_lowres = (0.01)^2*((Nr/Nr_lowres)^2)*Γprior_lowres[2:N_lowres+2,2:N_lowres+2];


# 
# W_stop = ones(Cdouble,N);
W_stop_lowres = ones(Cdouble,N_lowres);
τ0 = 1.0e1
# x00 = 0.5ones(Cdouble,N);
x00 = 0.5ones(Cdouble,N_lowres);
N_max_iter = 200000 # 0;
r_n_tol=0.0001;
r_y_tol=0.005;

# for each noise sample
# ρ_est_block    = zeros(Cdouble,Nnoise,Nr);
# ρ_est_cp_block = zeros(Cdouble,Nnoise,Nr);
ρ_est_block       = zeros(Cdouble,Nnoise,Nr_lowres);
ρ_est_cp_block    = zeros(Cdouble,Nnoise,Nr_lowres);
ρ_est_cp_block_un = zeros(Cdouble,Nnoise,Nr_lowres);
# ΓHΓyinv           = zeros(Cdouble,N_lowres,N_lowres);



for i in 1:Nnoise # does not look like this look is safe for multithreading, probably because of the shared memory Threads.@threads  
    println(i,"/",Nnoise)
    local Y = [y_tilde[i,:]; yd];

    local ρ_trunc = inv(Htrunc'*Htrunc)*Htrunc'*Y;

    # ρ_est_block[i,:] = [ρA_1[1]; ρ_trunc; ρA_1[N_lowres+2:end]]
    global ρ_est_block[i,:] = [ρA_1[1]; ρ_trunc; ρA_1[end]*ones(Cdouble,Nr_lowres-N0_lowres)]

    # local w = [(1.0/σnoise[i])*ones(Cdouble,Ndata); (1.0/0.1)*ones(Cdouble,N+1)];
    
    # local ρ_est,sn,taun,X_ALL,S_ALL,T_ALL,N_last = alg2_cp_gaussian(0.5ones(Cdouble,N),Y,Htrunc,w,W_stop;tau0=τ0,Niter=N_max_iter,r_n_tol=r_n_tol,r_y_tol=r_y_tol);
    # local ρ_est,sn,taun,X_ALL,S_ALL,T_ALL,N_last = alg2_cp_quad(0.5ones(Cdouble,N),y_tilde[i,:],yd,Htrunc,ΓI[:,:,i],Γd::Array{Cdouble,2},W_stop;τ0=τ0,Niter=N_max_iter,r_n_tol=r_n_tol,r_y_tol=r_y_tol)
    local ρ_est,sn,taun,X_ALL,S_ALL,T_ALL,N_last = alg2_cp_quad(0.5ones(Cdouble,N_lowres),y_tilde[i,:],yd,Htrunc,ΓI[:,:,i],Γd_lowres,W_stop_lowres;τ0=τ0,Niter=N_max_iter,r_n_tol=r_n_tol,r_y_tol=r_y_tol)
    println(N_last,"/",N_max_iter," i = ",i)

    # marginalization of uncertainty
    local ΓHΓyinv = zeros(Cdouble,N_lowres,N_lowres);
    for k in 1:Ndata
        # ΓHΓyinv = ΓHΓyinv + diagm(diag(ΓH[2:N0_lowres,2:N0_lowres,k])/ΓI[k,k,i])
        ΓHΓyinv = ΓHΓyinv + ΓH[2:N0_lowres,2:N0_lowres,k]/(ΓI[k,k,i]+(σB*σεH[k])^2)
    end
    local ρ_est_un,sn,taun,X_ALL,S_ALL,T_ALL,N_last = alg2_cp_quad_un(0.5ones(Cdouble,N_lowres),y_tildeμ[i,:],yd,Htrunc_un,ΓI[:,:,i],Γd_lowres,ΓHΓyinv,W_stop_lowres;τ0=τ0,Niter=N_max_iter,r_n_tol=r_n_tol,r_y_tol=r_y_tol)
    println(N_last,"/",N_max_iter," i = ",i)

    # ρ_est_cp_block[i,:] = [ρA_1[1]; ρ_est; ρA_1[N_lowres+2:end]]
    global ρ_est_cp_block[i,:] = [ρA_1[1]; ρ_est; ρA_1[end]*ones(Cdouble,Nr_lowres-N0_lowres)]
    global ρ_est_cp_block_un[i,:] = [ρA_1[1]; ρ_est_un; ρA_1[end]*ones(Cdouble,Nr_lowres-N0_lowres)]
end


μρ = dropdims(mean(ρ_est_cp_block,dims=1),dims=1);
Γρ = cov(ρ_est_cp_block);

μρ_un = dropdims(mean(ρ_est_cp_block_un,dims=1),dims=1);
Γρ_un = cov(ρ_est_cp_block_un);

figure()
plot(μ0.-r_lowres,ρ_est_block')
plot(μ0.-r,ρA_1)

figure()
# plot(ρ_est_block')
plot(μ0.-r_lowres,ρ_est_cp_block')
plot(μ0.-r,ρA_1)

figure()
# plot(ρ_est_block')
plot(μ0.-r_lowres,ρ_est_cp_block_un')
plot(μ0.-r,ρA_1)


figure()
plot(1000.0(r_lowres.-μ0),reverse(μρ),color="tab:blue",label="mean value")
fill_between(1000.0(r_lowres.-μ0),reverse(μρ-sqrt.(diag(Γρ))),reverse(μρ+sqrt.(diag(Γρ))),alpha=0.5,color="tab:blue",label="uncertainty")
plot(1000.0(r_lowres.-μ0),reverse(μρ_un),color="tab:orange",label="mean value marginal")
fill_between(1000.0(r_lowres.-μ0),reverse(μρ_un-sqrt.(diag(Γρ_un))),reverse(μρ_un+sqrt.(diag(Γρ_un))),alpha=0.5,color="tab:orange",label="uncertainty marginal")
plot(1000.0(r.-μ0),reverse(ρA_1),color="tab:green",label="GT")
legend(fontsize=14)
xlabel("depth [nm]",fontsize=14)
xticks(fontsize=14)
ylabel("concentration [a.u.]",fontsize=14)
yticks(fontsize=14)

# savefig("truncated_model_hard_values_50_datapoints.png")
# savefig("truncated_model_hard_values_50_datapoints.pdf")


# savefig("truncated_model_hard_values_50_datapoints_noise_3e-2.png")
# savefig("truncated_model_hard_values_50_datapoints_noise_3e-2.pdf")

# savefig("truncated_model_hard_values_20_datapoints.png")
# savefig("truncated_model_hard_values_20_datapoints.pdf")

# savefig("truncated_model_hard_values_10_datapoints.png")
# savefig("truncated_model_hard_values_10_datapoints.pdf")

# savefig("truncated_model_hard_values_10_datapoints_noise_1e-2.png")
# savefig("truncated_model_hard_values_10_datapoints_noise_1e-2.pdf")

# savefig("truncated_model_hard_values_10_datapoints_noise_1e-1.png")
# savefig("truncated_model_hard_values_10_datapoints_noise_1e-1.pdf")

# savefig("truncated_model_hard_values_5_datapoints.png")
# savefig("truncated_model_hard_values_5_datapoints.pdf")


include("postCovEst.jl")