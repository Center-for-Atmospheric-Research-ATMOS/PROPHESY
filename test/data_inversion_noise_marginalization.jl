
W_stop_lowres = ones(Cdouble,N_lowres);
τ0 = 1.0e1 # 
x00 = 0.5ones(Cdouble,N_lowres); # since the concentration is normalized by the bulk concentration, the initial state is taken as uniform with value 1/2
N_max_iter = 20000#00; #NOTE: for very low noise levels, the convergence can be fairly slow, so, one might consider increasing the maximum number of iterations
r_n_tol=0.001;
r_y_tol=0.001;
r_y_tol_un=r_y_tol; 

w = σw*ones(Cdouble,N_lowres); # not optimal because we know that the concentration varies more in the region near the surface rather than deep in the sample
Γsqrt = real(sqrt(corrCovariance(w;cor_len=10.0)));
Ns      = 1000000;
Ns_burn =  100000;


# for each noise sample
ρ_cp    = zeros(Cdouble,Nsample,Nr_lowres);
μρ_HI = zeros(Cdouble,N_lowres,Nsample);
Γρ_HI = zeros(Cdouble,N_lowres,N_lowres,Nsample);
deltaU = zeros(Cdouble,Ns,Nsample);

Threads.@threads  for i in 1:Nsample
    println(i,"/",Nsample)
    # argmax estimate
    local x00_loc = copy(x00);
    local yd_loc = copy(yd);
    local Htrunc_loc = copy(Htrunc);
    local ΓItrunc_loc = copy(ΓItrunc);
    local Γd_lowres_loc = copy(Γd_lowres);
    local W_stop_lowres_loc = copy(W_stop_lowres);
    local ρ_est,_,N_last = alg2_cp_quad_LM(x00_loc,y_tilde[i,:],yd_loc,Htrunc_loc,ΓItrunc_loc,Γd_lowres_loc,W_stop_lowres_loc;τ0=τ0,Niter=N_max_iter,r_n_tol=r_n_tol,r_y_tol=r_y_tol)
    global ρ_cp[i,:] = [ρA_1[1]; ρ_est; ρA_1[end]*ones(Cdouble,Nr_lowres-N0_lowres)];
    println(N_last,"/",N_max_iter)

    # posterior covariance estimation (conditional to data and model)
    println("posterior sampling: ",i)
    local Γsqrt_loc = copy(Γsqrt);
    local ΓIinv_loc = copy(ΓIinv);
    local Γd_lowres_inv_loc = copy(Γd_lowres_inv);
    local H_tilde_loc = copy(H_tilde);
    local D_tilde_loc = copy(D_tilde);
    μρ_HI[:,i],Γρ_HI[:,:,i],deltaU[:,i] = samplePosteriorMeanAndCov(ρ_cp[i,2:N0_lowres],Γsqrt_loc,y_tilde[i,:],yd_loc,ΓIinv_loc,Γd_lowres_inv_loc,H_tilde_loc,D_tilde_loc;Ns=Ns,Nburn=Ns_burn);
end

# variability due to noise in data
# marginalization P(ρ|H) = ∫P(ρ|H,y)P(y)dy -> mean and covariance of ρ_HI = argmax P(ρ|H,y)
mean_ρ_H = dropdims(mean(ρ_cp[1:Nsample,:],dims=1),dims=1);
var_ρ_H  = cov(ρ_cp[1:Nsample,:]);

# marginalization P(ρ|H) = ∫P(ρ|H,y)P(y)dy -> mean and covariance of ρ∼P(ρ|H)
μρ_H = dropdims(mean(μρ_HI,dims=2),dims=2); # μ_{ρ|H}
Γρ_H = dropdims(mean(Γρ_HI,dims=3),dims=3); # Γ_{ρ|H}




