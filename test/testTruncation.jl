## load the packages used in the estimation
# plotting
using PyPlot
fm = PyPlot.matplotlib.font_manager.json_load("/home/matthew/.cache/matplotlib/fontlist-v310.json") # TODO: look for the path automatically
# fm = PyPlot.matplotlib.font_manager.json_load("/home/mattoz/.cache/matplotlib/fontlist-v310.json")
# fm.findfont("serif", rebuild_if_missing=false)
# fm.findfont("serif", fontext="afm", rebuild_if_missing=false)
# rc("font",family="serif",serif="Computer Modern Roman")
rc("text", usetex=true)
rc("figure",max_open_warning=50)
using myPlot

# data manipulation (loading, writing, etc)
using Printf

# scientific package from the official Julia repositories
using LinearAlgebra
# using Statistics
# using DSP
# using SpecialMatrices
# using Polynomials
using StatsBase

# implemented scientific packages
using utilsFun  # for the softMax functions

# modeling XPS
using XPSpack
using XPSinv



## acquisition setup
ħν = 900.0;
μKe = ħν-285.0;
α = 1.0
T = 1.0
Fν = 1.0;
Nke = 200;
Ke = collect(range(μKe-2.0,μKe+2.0,length=Nke));
dKe = Ke[2]-Ke[1]
σν0 = 0.6;
σν = σν0*((0.7/sqrt(2π*0.2^2))*exp.(-0.5*(Ke.-(μKe-0.5)).^2/0.2^2) .+ (0.3/sqrt(2π*0.5^2))*exp.(-0.5*(Ke.-(μKe+0.5)).^2/0.5^2));
λe0 = 0.002;

wsAcq = XPSacq(ħν,μKe,α,T,Fν,Ke,σν,λe0);

## geometry setup
k0 = 10; # 5
Nr = 101;  # 51;
Nθ = 256;
Ny = 256;
μ0 = 20.0 # 100.0;
L = 200.0;
x0 = sqrt(2.0)*100.0
y0 = 0.0;
z0 = 100.0
r = collect(range(μ0-k0*λe0,μ0,length=Nr))
θ0 = atan(x0,z0)
θ = collect(range(θ0-π/2.0,θ0+π/2.0,Nθ));
y = collect(range(-L/2.0,L/2.0,length=Ny));

wsGeom = cylinderGeom(x0,y0,z0,μ0,r,θ,y);



# simulate some data (one point in the kinetic energy spectrum for four different concentration profiles)
ρA_1 = logistic.(1000.0reverse(μ0.-r).-2.0,0.0,1.0,2.0);
# ρA_1 = logistic.(1000.0reverse(μ0.-r).-2.0,0.0,1.0,2.0) .+ 2.0exp.(-(1000.0reverse(μ0.-r).-1.0).^2. /(2.0*0.25^2));
# ρA_1 = logistic.(1000.0reverse(μ0.-r).-2.0,0.0,1.0,2.0) .+ exp.(-(1000.0reverse(μ0.-r).-1.5).^2. /(2.0*0.5^2));
# ρA_1 = exp.(-(1000.0reverse(μ0.-r).-2.5).^2. /(2.0*0.5^2));

# covariance matrix for the a priori distribution



Γprior = zeros(Cdouble,Nr,Nr);
cor_len = 5.0;
for i in 1:Nr
    # Γprior[i,i] = (1.0-ρA_1[i]+0.1)^2;
    Γprior[i,i] =  1.0
    for j in i+1:Nr
        Γprior[i,j] = Γprior[i,i]*exp(-(i-j)^2/(0.5*cor_len^2))
        Γprior[j,i] = Γprior[i,j]
    end
end

# figure(); imshow(Γprior); colorbar()


Dprior = D2nd(Nr) # D2nd(Nr+2)[:,2:end-1];
# Bprior = 1.0e-12Dprior'*inv(Γprior[2:end-1,2:end-1])*Dprior; # 1.0e-8
Bprior = 1.0e-12Dprior'*inv(Γprior[2:end-1,2:end-1])*Dprior;


# measurement operator (only the geometical term since the other comes as a multiplicative scalar estimated from the data)
Ndata = 6 # 25
H_better = zeros(Cdouble,Ndata,Nr);
# λbetter0  = 1.0e-3*[1.0; 1.5; 2.0; 2.5; 3.0]; # these are some eal values that would nice to be able to access... but that would not be sufficient to make the uncertainty small enough
λbetter0 = 1.0e-3collect(range(1.3,2.5,Ndata));

for i in 1:Ndata
    H_better[i,:],_,_,_,_ = cylinder_gain_H(r,θ,y,x0,y0,z0,μ0,λbetter0[i]);
end

H_better = reverse(H_better,dims=2); #
figure(); plot(r.-μ0,H_better'); plot(r.-μ0,ρA_1)
plot(r[end-(29+50):end].-μ0,ρA_1[end-(29+50):end])


# generate some data (data point and covariance)
Nnoise = 5;
σnoise = [0.001; 0.01; 0.1; 1.0; 10.0];
Nnoise = 10;
σnoise = 0.001*ones(Cdouble,Nnoise);

y_data = zeros(Cdouble,Nnoise,Ndata);
ΓI = zeros(Cdouble,Ndata,Ndata,Nnoise);
ΓIsqrt = zeros(Cdouble,Ndata,Ndata,Nnoise);
detΓI = zeros(Cdouble,Nnoise);
ΓIinv = zeros(Cdouble,Ndata,Ndata,Nnoise);
for i in 1:Nnoise
    ΓI[:,:,i] = σnoise[i]^2*diagm(ones(Cdouble,Ndata)); # iid data noise
    ΓIsqrt[:,:,i] = sqrt(ΓI[:,:,i]);
    detΓI[i] = det(ΓI[:,:,i]);
    ΓIinv[:,:,i] = inv(ΓI[:,:,i]);
    y_data[i,:] = H_better*ρA_1 + ΓIsqrt[:,:,i]*randn(Cdouble,Ndata);
end
y_data[y_data.<0.0] = -y_data[y_data.<0.0];

Δy = H_better[:,end-(29+50):end]*ρA_1[end-(29+50):end];
δy = H_better[:,1]*ρA_1[1];

figure(); plot(y_data'); plot((y_data.-Δy')'); plot((y_data.-(Δy+δy)')')
plot(Δy); plot(δy)

y_tilde = y_data.-(Δy+δy)';
H_tilde = H_better[:,2:end-(29+50)]

N = size(H_tilde,2);

D_2nd = D2nd(N);
# B

A = [H_tilde; Matrix{Cdouble}(I,N,N); 10000.0D_2nd]; #  B];
w = (1.0/σnoise[1])*ones(Cdouble,Ndata);
W_stop = ones(Cdouble,N);
τ0 = 1.0e1 # 4
x00 = 0.5ones(Cdouble,N); # 0.0*ρA_1[idx_res]
s0 = A*x00;
N_max_iter = 2000#0 # 00;
r_n_tol=1.0
r_y_tol=0.005;


# xn,sn,taun,N_last = alg2_cp_gaussian_un_no_mem(x0,s0,IA_1,A,γ,γ_H,γ_D,W_stop;tau0=τ0,Niter=N_max_iter,r_n_tol=r_n_tol,r_y_tol=r_y_tol);
# xn,sn,taun,F_ALL,G_ALL,F_STAR_ALL,G_STAR_ALL,INNER_PROD,T_ALL,N_last= alg2_cp_gaussian_un_no_mem(x0,s0,IA_1,A,γ,γ_H,γ_D,W_stop;tau0=τ0,Niter=N_max_iter,r_n_tol=r_n_tol,r_y_tol=r_y_tol);

# Ns = 5 # 00;
# ρ_samples = zeros(Cdouble,Ns,N);

# ρ_samples[1,:],_,_,_,_,_,_,_,_,N_last= alg2_cp_gaussian_un_no_mem_val(x0,s0,IA_1_clean,ρ0,ρB,A,γ,γ_H,γ_D,γ0,γB,W_stop;tau0=τ0,Niter=N_max_iter,r_n_tol=r_n_tol,r_y_tol=r_y_tol);

ρ_est,sn,taun,X_ALL,S_ALL,T_ALL,N_last = alg2_cp_gaussian(x00,y_tilde[1,:],A,w,W_stop;tau0=τ0,Niter=N_max_iter,r_n_tol=r_n_tol,r_y_tol=r_y_tol)

figure(); plot([ρA_1[1]; ρ_est; ρA_1[end-(29+50):end]]); plot(ρA_1)