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
k0 = 5;
Nr = 51;
Nθ = 256;
Ny = 256;
μ0 = 100.0;
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


# measurement operator (only the geometical term since the other comes as a multiplicative scalar estimated from the data)
Ndata = 6 # 25
H_better = zeros(Cdouble,Ndata,Nr);
# λbetter0  = 1.0e-3*[1.0; 1.5; 2.0; 2.5; 3.0]; # these are some eal values that would nice to be able to access... but that would not be sufficient to make the uncertainty small enough
λbetter0 = 1.0e-3collect(range(1.3,2.5,Ndata));

for i in 1:Ndata
    H_better[i,:],_,_,_,_ = cylinder_gain_H(r,θ,y,x0,y0,z0,μ0,λbetter0[i]);
end

H_better = reverse(H_better,dims=2); #
figure(); plot(r.-μ0,H_better')

# generate some data (data point and covariance)
# Nnoise = 5;
# σnoise = [0.001; 0.01; 0.1; 1.0; 10.0];
Nnoise = 10;
σnoise = 1.0*ones(Cdouble,Nnoise);
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



figure(); plot(y_data')

# 0.001 -> 1.0e8
# 0.01  -> 1.0e6
# 0.1   -> 1.0e4
# 1.0   -> 1.0e2
# 10.0  -> 1.0e0





##
## the sampling process: sample the a posteriori distribution around a maximum to estimate a covariance matrix which represents the uncertainty to be expected when using this Bayesian model
##

# square root matrix of the generative covariance matrix (the covariance in the distribution used for generating new samples)
σw = 1.0e-1 # 0.001; # small compared with the amplitude of the state 
w = σw*ones(Cdouble,Nr); # not optimal because we know that the concentration varies more in the region near the surface rather than deep in the sample
Γsqrt = real(sqrt(corrCovariance(w;cor_len=15.0)));
p0 = 0.05 # 0.02; # starting acceptance rate of uphill moves
# ρB = [ρA_1[1]; ρA_1[end]]; # known values 
# σB = [0.01; 0.01];         # how much do we trust these known values
# for noise levels: σnoise = [0.001; 0.01; 0.1; 1.0; 10.0], use the following entropy weights
# wE = [1.0e8; 1.0e6; 1.0e4; 1.0e2; 1.0e0]
wE = 1.0e2*ones(Cdouble,Nnoise);
ρE = logistic.(1000.0reverse(μ0.-r).-2.0,0.0,1.0,2.0); # purposefully choose a wrong a priori profile
Ns = 1000000; # number of samples... no idea a priori how many samples are needed

Γρ_IG = zeros(Cdouble,Nr,Nr,Nnoise);
μρ_IG = zeros(Cdouble,Nr,Nnoise);

PLOT_FIG = false

for k in 1:Nnoise
    # actually run the sampling
    # ρ_all = samplePosterior(zeros(Cdouble,Nr),Γsqrt,p0*ones(Cdouble,Ns),y_data[k,:],ΓIinv[:,:,k],H_better,Bprior,ρB,σB;Ns=Ns,psmooth=1.999);
    ρ_all = samplePosteriorEntropy(zeros(Cdouble,Nr),Γsqrt,p0*ones(Cdouble,Ns),y_data[k,:],ΓIinv[:,:,k],H_better,ρE,wE[k];Ns=Ns)

    # compute a covariance matrix from the samples 
    μρ_IG[:,k] = dropdims(mean(ρ_all,dims=1),dims=1);
    Γρ_IG[:,:,k] = cov(ρ_all);

    if PLOT_FIG
        figure();
        plot(μ0.-r,ρA_1,label="GT")
        plot(μ0.-r,μρ_IG[:,k],color="blue",label="average")
        fill_between(μ0.-r,ρA_1-sqrt.(diag(Γρ_IG[:,:,k])),ρA_1+sqrt.(diag(Γρ_IG[:,:,k])),alpha=0.5,color="tab:blue",label="uncertainty")
        # figure(); imshow(Γρ_IG[:,:,k]); colorbar()
        # figure(); plot(μ0.-r,sqrt.(diag(Γρ_IG[:,:,k])))


        
        # observe the quantities used for the sampling (note that the burn in period is not too long with the communication kernel in use)
        Elikelihood  = zeros(Cdouble,Ns+1);
        EpriorEntropy = zeros(Cdouble,Ns+1);
        for i in 1:Ns+1
            Elikelihood[i]   = (y_data[k,:]-H_better*ρ_all[i,:])'*ΓIinv[:,:,k]*(y_data[k,:]-H_better*ρ_all[i,:])
            EpriorEntropy[i] = sum(ρ_all[i,:]-ρE - ρ_all[i,:].*log.(ρ_all[i,:]./ρE));
        end
        Etot = Elikelihood+EpriorEntropy;
        Etot[isnan.(Etot)] .= Inf;

        figure()
        semilogx(collect(1:Ns+1),Elikelihood,label="likelihood")
        semilogx(collect(1:Ns+1),EpriorEntropy,label="entropy a priori")
        semilogx(collect(1:Ns+1),Etot,label="total")
        legend()

        val_min,idx_min = findmin(Etot)
        figure();
        plot(μ0.-r,ρA_1,label="GT")
        plot(μ0.-r,ρ_all[idx_min,:],label="min Etot")
        plot(μ0.-r,μρ_IG[:,k],color="blue",label="average")
        fill_between(μ0.-r,ρA_1-sqrt.(diag(Γρ_IG[:,:,k])),ρA_1+sqrt.(diag(Γρ_IG[:,:,k])),alpha=0.5,color="tab:blue",label="uncertainty")
    end
    
end


figure()
plot(μ0.-r,ρA_1,label="GT")
μ_mean = dropdims(median(μρ_IG,dims=2),dims=2);
Γ_mean = dropdims(median(Γρ_IG,dims=3),dims=3);
plot(μ0.-r,μ_mean,color="blue",label="average")
fill_between(μ0.-r,ρA_1-sqrt.(diag(Γ_mean)),ρA_1+sqrt.(diag(Γ_mean)),alpha=0.5,color="tab:blue",label="uncertainty")
legend(fontsize=14)
xlabel("distance [\$\\mu\$m]",fontsize=14)
ylabel("relative concentration [a.u.]",fontsize=14)

# savefig("rho1_posterior_cov_noise_1e-2_entropy.png")
# savefig("rho1_posterior_cov_noise_1e-2_entropy.pdf")

# savefig("rho1_posterior_cov_noise_1e-1_entropy.png")
# savefig("rho1_posterior_cov_noise_1e-1_entropy.pdf")

# savefig("rho1_posterior_cov_noise_1e0_entropy.png")
# savefig("rho1_posterior_cov_noise_1e0_entropy.pdf")


[norm(Γρ_IG[:,:,i]) for i in 1:Nnoise]



# Eval  = zeros(Cdouble,Ns);
# Eval2 = zeros(Cdouble,Ns);
# for i in 1:Ns
#     new_state = ρA_1 + Γsqrt*randn(Cdouble,Nr);
#     new_state[new_state.<=0.0] .= 1.0e-6
#     Eval2[i]   = (y_data[1,:]-H_better*new_state)'*ΓIinv[:,:,1]*(y_data[1,:]-H_better*new_state)
#     Eval[i] = sum(new_state-ρE) - sum(new_state.*log.(new_state./ρE))
# end
# figure(); hist(Eval,100)
# figure(); hist(Eval2,100)

# figure(); plot(ρA_1 + Γsqrt*randn(Cdouble,Nr))

# Eval[isnan.(Eval)] .= Inf;

# if false
#     # explaining the difference between entropy regularization and smoothness regularization:
#     x1 = collect(range(0.1,5.0,length=50));
#     x2 = collect(range(0.1,5.0,length=50));
#     x1E = 1.5; x2E = 0.5;
#     Γxx = inv(2.0*[1.0 0.25; 0.25 0.1]);

#     y1 = 3.5; y2 = 3.5;
#     σ1 = 1.0; σ2 = 1.0;
#     Plikelihood = exp.(-0.5*((y1.-x1)./σ1).^2)*exp.(-0.5*((y2.-x2)./σ2).^2)';

#     figure()
#     imshow(Plikelihood)
#     colorbar()

#     Psmooth  = zeros(Cdouble,50,50);
#     Pentropy = zeros(Cdouble,50,50);
#     for i1 in 1:50
#         for i2 in 1:50
#             xx = [x1[i1]; x2[i2]];
#             Psmooth[i1,i2] = exp(-0.5xx'Γxx*xx)
#             Pentropy[i1,i2] = exp(2.0e0*(x1[i1] - x1E - x1[i1]*log(x1[i1]/x1E) + x2[i2] - x2E - x2[i2]*log(x2[i2]/x2E)))
#         end
#     end

#     figure(figsize=[10,5])
#     ax1 = subplot(121)
#     contour(x1,x2,Plikelihood,label="likelihood")
#     contour(x1,x2,Psmooth,label="smoothness a priori")
#     contour(x1,x2,Plikelihood.*Psmooth,label="posterior")
#     # title("Smoothness")
#     xlabel("x1")
#     ylabel("x2")

#     ax2 = subplot(122)
#     contour(x1,x2,Plikelihood)
#     contour(x1,x2,Pentropy)
#     contour(x1,x2,Plikelihood.*Pentropy)
#     # title("Entropy")
#     xlabel("x1")
#     ylabel("x2")

#     tight_layout(pad=1.0, w_pad=0.5, h_pad=0.2)


#     # ax1.annotate("a)", xy=(3, 1),  xycoords="data", xytext=(-1.27, 0.99), textcoords="axes fraction", color="black",fontsize=14)
#     ax1.annotate("a)", xy=(3, 1),  xycoords="data", xytext=(-0.08, 0.99), textcoords="axes fraction", color="black",fontsize=14)
#     ax2.annotate("b)", xy=(3, 1),  xycoords="data", xytext=(-0.08, 0.99), textcoords="axes fraction", color="black",fontsize=14)

#     ax1.annotate("Likelihood", xy=(3, 1),  xycoords="data", xytext=(0.6, 0.9), textcoords="axes fraction", color="black",fontsize=14)
#     ax1.annotate("a priori", xy=(3, 1),  xycoords="data", xytext=(0.01, 0.05), textcoords="axes fraction", color="black",fontsize=14)
#     ax1.annotate("posterior", xy=(3, 1),  xycoords="data", xytext=(0.08, 0.55), textcoords="axes fraction", color="black",fontsize=14)
#     ax1.annotate("Smoothness", xy=(3, 1),  xycoords="data", xytext=(0.6, 0.1), textcoords="axes fraction", color="black",fontsize=14)

#     ax2.annotate("Likelihood", xy=(3, 1),  xycoords="data", xytext=(0.6, 0.9), textcoords="axes fraction", color="black",fontsize=14)
#     ax2.annotate("a priori", xy=(3, 1),  xycoords="data", xytext=(0.03, 0.17), textcoords="axes fraction", color="black",fontsize=14)
#     ax2.annotate("posterior", xy=(3, 1),  xycoords="data", xytext=(0.25, 0.45), textcoords="axes fraction", color="black",fontsize=14)
#     ax2.annotate("Entropy", xy=(3, 1),  xycoords="data", xytext=(0.6, 0.1), textcoords="axes fraction", color="black",fontsize=14)

#     # savefig("smoothness_vs_entropy_formulation.png")
#     # savefig("smoothness_vs_entropy_formulation.pdf")




#     figure()
#     imshow(Psmooth.*Plikelihood)
#     colorbar()


#     figure()
#     imshow(Pentropy.*Plikelihood)
#     colorbar()


# end