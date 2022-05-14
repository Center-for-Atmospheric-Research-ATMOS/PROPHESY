ρ_cp_HI        = zeros(Cdouble,N_model_sample,Nr_lowres);
deltaUh        = zeros(Cdouble,Ns,N_model_sample);
μρ_HI_sample   = zeros(Cdouble,N_lowres,N_model_sample);
Γρ_HI_sample   = zeros(Cdouble,N_lowres,N_lowres,N_model_sample);
# 
Threads.@threads for m in 1:N_model_sample # this loop is just meant for sampling the model, so, each iteration is independent
    println(m,"/",N_model_sample)
    local num_sample;
    if (m<10)
        num_sample = string("00",m,"/")
    elseif ((m>=10) & (m<100))
        num_sample = string("0",m,"/")
    else
        num_sample = string(m,"/")
    end
    local dfH_lowres = CSV.File(string(model_folder_lowres_un,"/",num_sample,"H_lowres.csv");header=true,ntasks=1) |> DataFrame
    local H_lowres   = Matrix{Cdouble}(dfH_lowres);



    # slice the model (3 terms: boundary, surface and bulk)
    local H0 = H_lowres[:,1];
    local H_tilde = H_lowres[:,2:N0_lowres];
    local Hb = H_lowres[:,N0_lowres+1:end];

    # data and noise correction 
    local Δy      = dropdims(sum(Hb,dims=2),dims=2)*ρA_1[end];
    local δy      = H0*ρA_1[1];
    local y_tilde = repData.-(Δy+δy)';
    local Hnot    = [H0 Hb];
    local Hnot1   = sum(Hnot;dims=2);
    local ΓItrunc = ΓI + σB^2*Hnot1*Hnot1';
    local ΓIinv   = inv(ΓItrunc);


    # smoosh together the several part of the model into augmented operators
    local Htrunc = [H_tilde; D_tilde];                                     # conditional to data and measurement model

    # data inversion: estimation of the concentration profile using CP algorithm
    local x00_loc = copy(x00);
    local yd_loc = copy(yd);
    local Γd_lowres_loc = copy(Γd_lowres);
    local W_stop_lowres_loc = copy(W_stop_lowres);
    println(i_sample," data, model: ",m)
    local ρ_est,_,N_last = alg2_cp_quad_LM(x00_loc,y_tilde[i_sample,:],yd_loc,Htrunc,ΓItrunc,Γd_lowres_loc,W_stop_lowres_loc;τ0=τ0,Niter=N_max_iter,r_n_tol=r_n_tol,r_y_tol=r_y_tol)
    global ρ_cp_HI[m,:] = [ρA_1[1]; ρ_est; ρA_1[end]*ones(Cdouble,Nr_lowres-N0_lowres)]
    println(N_last,"/",N_max_iter," model ",m)

    # sample the posterior distribution for the fixed data y_tilde[i_sample,:]
    # posterior covariance estimation (conditional to data and model)
    local Γsqrt_loc = copy(Γsqrt);
    local Γd_lowres_inv_loc = copy(Γd_lowres_inv);
    local H_tilde_loc = copy(H_tilde);
    local D_tilde_loc = copy(D_tilde);
    println("posterior sampling: ",i_sample," model ",m)
    μρ_HI_sample[:,m],Γρ_HI_sample[:,:,m],deltaUh[:,m] = samplePosteriorMeanAndCov(ρ_cp_HI[m,2:N0_lowres],Γsqrt_loc,y_tilde[i_sample,:],yd_loc,ΓIinv,Γd_lowres_inv_loc,H_tilde_loc,D_tilde_loc;Ns=Ns,Nburn=Ns_burn);
end

# compute the estimate variability due to the variation in the model
mean_ρ_y = dropdims(mean(ρ_cp_HI[1:N_model_sample,:],dims=1),dims=1);
var_ρ_y = cov(ρ_cp_HI[1:N_model_sample,:]);

# marginal mean and covariance
μρ_y = dropdims(mean(μρ_HI_sample,dims=2),dims=2); # μ_{ρ|y}


Γρ_y_mean = dropdims(mean(Γρ_HI_sample,dims=3),dims=3)
ΔΓρ_y1 = (1.0/N_model_sample^2)*μρ_y*μρ_y';
ΔΓρ_y2 = μρ_HI_sample[:,1]*μρ_HI_sample[:,1]';
[global ΔΓρ_y2 = ΔΓρ_y2 + μρ_HI_sample[:,i]*μρ_HI_sample[:,i]'; for i in 2:N_model_sample]
ΔΓρ_y2 = (1.0/N_model_sample)*ΔΓρ_y2;

Γρ_y = Γρ_y_mean - ΔΓρ_y1 + ΔΓρ_y2; # Γ_{ρ|y}

