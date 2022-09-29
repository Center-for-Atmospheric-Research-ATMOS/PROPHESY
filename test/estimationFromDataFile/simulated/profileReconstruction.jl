##
## inversion model
##

# standard deviation of the known values
σB     = 0.1;        

# standard deviation of the smoothness a priori
σd     = 0.05 # 0.1;
if FLAG_0001
    σd     = 0.025
end
if FLAG_0002
    σd     = 0.05
end
if FLAG_0003
    σd     = 0.075
end
if FLAG_0004
    σd     = 0.2
end
cor_len_lowres = 2.5;

# amplitude of the communication mecanism for sampling the a posteriori model
σw     = 5.0e-4;  # only used for sampling
if FLAG_0001
    σw     = 5.0e-4
end
if FLAG_0002
    σw     = 5.0e-4
end
if FLAG_0003
    σw     = 1.0e-3
end
if FLAG_0004
    σw     = 1.0e-3
end

# deeper than some distance, the signal is not likely to be disantangled
d0 = 2.5e-3 # [μm] maximum depth
N0 = findlast(r.-μ0.<=-d0);
N = Nr-N0;

# slice the model (3 terms: boundary, surface and bulk)
H0 = Hgeom[:,end];
H_tilde = Hgeom[:,N0:end-1];
Hb = Hgeom[:,1:N0-1];
Hnot = [Hb H0];
Hnot1 = sum(Hnot;dims=2);


# data correction
Δy = dropdims(sum(Hb,dims=2),dims=2); # in the bulk (deeper than 5 nm) if the data are normalized by the bulk concentration, then the normalized concentration in the bulk is 1, otherwise ρB
if FLAG_0004
    Δy = 0.0*Δy
end
δy = H0*0.0;                          # outside the sample the concentration is nearly 0
y_tilde  = y_data_1-(Δy+δy);          # total correction



# regularization (smoothness: applied as sparsity in the second order difference)
DN = D2nd(N+3);
D0 = DN[:,end];
D_tilde = DN[:,3:end-1];
Db = DN[:,1:2];

# correction of the regularization "data"
Δyd = -dropdims(sum(Db,dims=2),dims=2); # if data not normalized by bulk concentration, multiply by ρB
if FLAG_0004
    Δyd = 0.0*Δyd
end
δyd = -D0*0.0;
yd = Δyd+δyd;


# smoosh together the several part of the model into augmented operators
Htrunc = [H_tilde; D_tilde];                                     # conditional to data and measurement model

#
# covariances: the crafed Bayesian models assumes that some covariance are known, i.e. measurement noise, smoothness, known values and measurement operator (the last only in the marginalized case)
#

# measurement noise covariance
ΓI = diagm(σnoise.^2);
ΓItrunc = ΓI + σB^2*Hnot1*Hnot1';
ΓIinv = inv(ΓItrunc);

# covariance matrix for the a priori distribution (second order difference)
Γprior = zeros(Cdouble,Nr,Nr);
for i in 1:Nr
    Γprior[i,i] =  1.0
    for j in i+1:Nr
        Γprior[i,j] = Γprior[i,i]*exp(-(i-j)^2/(0.5*cor_len_lowres^2));
        Γprior[j,i] = Γprior[i,j];
    end
end
Γd = (N/NdataC1s)*(σd^2)*Γprior[N0-1:end-1,N0-1:end-1];  # scale the a priori strength with the quantity of data, so that it is possible to compare the results
Γd_inv = inv(Γd);


##
## reconstruction
##
W_stop = ones(Cdouble,N);
W_stop = collect(LinRange(1.0,10.0,N));
τ0 = 1.0e1 # 
x00 = 0.5ones(Cdouble,N); # since the concentration is normalized by the bulk concentration, the initial state is taken as uniform with value 1/2
N_max_iter = 200000# 0#00; 
r_n_tol=0.001;
r_y_tol=0.001;
r_y_tol_un=r_y_tol; 
ρ_cp    = zeros(Cdouble,Nr);
μρ_HI = zeros(Cdouble,N);
Γρ_HI = zeros(Cdouble,N,N);

ρ_est,_,N_last = alg2_cp_quad_LM(x00,y_tilde,yd,Htrunc,ΓItrunc,Γd,W_stop;τ0=τ0,Niter=N_max_iter,r_n_tol=r_n_tol,r_y_tol=r_y_tol);

ρ_cp = [ones(Cdouble,N0-1); ρ_est; 0.0];
if FLAG_0004
    ρ_cp = [zeros(Cdouble,N0-1); ρ_est; 0.0];
end


##
## sampling the a posteriori model
##
if SAMPLING
    w = σw*ones(Cdouble,N); # not optimal because we know that the concentration varies more in the region near the surface rather than deep in the sample
    Γsqrt = real(sqrt(corrCovariance(w;cor_len=10.0))); # 5.0
    Ns      = 1000000;
    Ns_burn =  100000;
    deltaU = zeros(Cdouble,Ns);

    global μρ_HI,Γρ_HI,deltaU = samplePosteriorMeanAndCov(ρ_cp[N0:end-1],Γsqrt,y_tilde,yd,ΓIinv,Γd_inv,H_tilde,D_tilde;Ns=Ns,Nburn=Ns_burn);

    global stdρ_HI = sqrt.(diag(Γρ_HI));
end


