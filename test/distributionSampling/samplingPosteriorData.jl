## load the packages used in the estimation
# plotting
using PyPlot
rc("text", usetex=true)
rc("figure",max_open_warning=50)

# data manipulation (loading, writing, etc)
using Printf

# scientific package from the official Julia repositories
using LinearAlgebra
using StatsBase


# modeling XPS
using XPSpack
using XPSsampling



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
# ρA_2 = logistic.(1000.0reverse(μ0.-r).-2.0,0.0,1.0,2.0) .+ 2.0exp.(-(1000.0reverse(μ0.-r).-1.0).^2. /(2.0*0.25^2));
# ρA_3 = logistic.(1000.0reverse(μ0.-r).-2.0,0.0,1.0,2.0) .+ exp.(-(1000.0reverse(μ0.-r).-1.5).^2. /(2.0*0.5^2));
# ρA_4 = exp.(-(1000.0reverse(μ0.-r).-2.5).^2. /(2.0*0.5^2));

# covariance matrix for the a priori distribution



Γprior = zeros(Cdouble,Nr,Nr)
cor_len = 5.0;
for i in 1:Nr
    # Γprior[i,i] = (0.005*(1.0-ρA_1[i]+0.1))^2;
    Γprior[i,i] =  0.005^2 
    for j in i+1:Nr
        Γprior[i,j] = Γprior[i,i]*exp(-(i-j)^2/(0.5*cor_len^2))
        Γprior[j,i] = Γprior[i,j]
    end
end

figure(); imshow(Γprior); colorbar()


Dprior = diagm(Nr-2,Nr,1 => 2ones(Cdouble,Nr-2), 0 => -ones(Cdouble,Nr-2) ,2 => -ones(Cdouble,Nr-2))
# Bprior = 1.0e8Dprior'*inv(Γprior[2:end-1,2:end-1])*Dprior; # 1.0e-8
Bprior = 1.0e-16Dprior'*inv(Γprior[2:end-1,2:end-1])*Dprior; # 1.0e-8
Cprior = inv(Bprior);
Dsqrt = real(sqrt(Cprior));




# measurement operator (only the geometical term since the other comes as a multiplicative scalar estimated from the data)
Ndata = 5
H_better = zeros(Cdouble,Ndata,Nr);
Γbetter  = zeros(Cdouble,Ndata,Nr,Nr);
λbetter0  = 1.0e-3*[1.0; 1.5; 2.0; 2.5; 3.0];
Nλ = 21;
λbetter = zeros(Cdouble,Ndata,Nλ);

for i in 1:Ndata
    λbetter[i,:] = collect(range(0.9995λbetter0[i],1.0005λbetter0[i],length=Nλ));
    Γbetter[i,:,:],H_better[i,:] = cov_H_cylinder(r,θ,y,x0,y0,z0,μ0,λbetter[i,:],(1.0/(λbetter[i,end]-λbetter[i,1]))*ones(Cdouble,Nλ));
end

H_better = reverse(H_better,dims=2);
figure(); plot(r.-μ0,H_better')


Γbetter = reverse(Γbetter,dims=(2,3));

for i in 1:Ndata
    figure(); imshow(Γbetter[i,:,:]); colorbar()
end

# generate some data (data point and covariance)
Nnoise = 6;
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


# square root matrix of the generative covariance matrix (the covariance in the distribution used for generating new samples)
σw = 10.e-2 # 0.001; # small compared with the amplitude of the state 
w = σw*ones(Cdouble,Nr); # not optimal because we know that the concentration varies more in the region near the surface rather than deep in the sample
# w = σw*(1.0.-ρA_1.+0.1); # too optimal because we known the solution (but the general sigmoid shape could be used because it's not a big a priori)
Γsqrt = real(sqrt(corrCovariance(w;cor_len=15.0)));
ρB = [ρA_1[1]; ρA_1[end]];
σB = [0.01; 0.01];
Ns = 100 # 1000000; #NOTE: 100 for quick test and 1000000 for estimation (but it may take long)

Γρ_I = zeros(Cdouble,Nr,Nr,Nnoise);
μρ_I = zeros(Cdouble,Nr,Nnoise);

PLOT_FIG = (false & (Threads.nthreads()==1)) # to much mess if more than one thread, generating a segfault

Threads.@threads  for k in 1:Nnoise
    println(k,"/",Nnoise)

    local ρ_all = samplePosteriorModelMargin(0.0ρA_1,Γsqrt,y_data[k,:],ΓIinv[:,:,k],H_better,Γbetter,Bprior,ρB,σB;Ns=Ns,psmooth=1.999);

    # compute a covariance matrix from the samples 
    global μρ_I[:,k] = dropdims(mean(ρ_all,dims=1),dims=1);
    global Γρ_I[:,:,k] = cov(ρ_all);

    if PLOT_FIG
        figure();
        plot(μ0.-r,ρA_1)
        plot(μ0.-r,μρ_I[:,k],color="blue")
        fill_between(μ0.-r,ρA_1-sqrt.(diag(Γρ_I[:,:,k])),ρA_1+sqrt.(diag(Γρ_I[:,:,k])),alpha=0.5,color="tab:blue",label="uncertainty")
        figure(); imshow(Γρ_I[:,:,k]); colorbar()
        figure(); plot(μ0.-r,sqrt.(diag(Γρ_I[:,:,k])))


        # TODO: observe the burn in period and don't use it for the computation of the mean and the covariance... not really showing up (which is good news)
        local Elikelihood  = zeros(Cdouble,Ns+1);
        local EpriorSmooth = zeros(Cdouble,Ns+1);
        local EpriorVal    = zeros(Cdouble,Ns+1);
        local EpriorModel  = zeros(Cdouble,Ns+1);
        for i in 1:Ns+1
            Elikelihood[i]  = (y_data[k,:]-H_better*ρ_all[i,:])'*ΓIinv[:,:,k]*(y_data[k,:]-H_better*ρ_all[i,:])
            EpriorSmooth[i] = ρ_all[i,:]'Bprior*ρ_all[i,:]
            EpriorVal[i] = ((ρ_all[i,1]  -ρB[1])^2)/(σB[1]^2) + ((ρ_all[i,end]  -ρB[2])^2)/(σB[2]^2)
            for g in 1:Ndata
                EpriorModel[i] = EpriorModel[i] + ΓIinv[g,g,k]*(ρ_all[i,:]'*Γbetter[g,:,:]*ρ_all[i,:])
            end
        end

        local Etot = Elikelihood+EpriorSmooth+EpriorVal+EpriorModel;
        Etot[isnan.(Etot)] .= Inf;
        local val_min,idx_min = findmin(Etot);

        figure()
        semilogx(collect(1:Ns+1),Elikelihood,label="likelihood")
        semilogx(collect(1:Ns+1),EpriorSmooth,label="smoothness a priori")
        semilogx(collect(1:Ns+1),EpriorVal,label="values a priori")
        semilogx(collect(1:Ns+1),EpriorModel,label="operator a priori")
        semilogx(collect(1:Ns+1),Etot,label="total")
        legend()

        figure();
        plot(μ0.-r,ρA_1,label="GT")
        plot(μ0.-r,ρ_all[idx_min,:],label="min Etot")
        plot(μ0.-r,μρ_I[:,k],color="blue",label="average")
        fill_between(μ0.-r,ρA_1-sqrt.(diag(Γρ_I[:,:,k])),ρA_1+sqrt.(diag(Γρ_I[:,:,k])),alpha=0.5,color="tab:blue",label="uncertainty")
    end

end


figure()
plot(μ0.-r,ρA_1,label="GT")
μ_mean = dropdims(mean(μρ_I,dims=2),dims=2)
Γ_mean = dropdims(mean(Γρ_I,dims=3),dims=3)
plot(μ0.-r,μ_mean,color="blue",label="average")
fill_between(μ0.-r,ρA_1-sqrt.(diag(Γ_mean)),ρA_1+sqrt.(diag(Γ_mean)),alpha=0.5,color="tab:blue",label="uncertainty")
legend(fontsize=14)
xlabel("distance [\$\\mu\$m]",fontsize=14)
ylabel("relative concentration [a.u.]",fontsize=14)

# savefig("rho1_posterior_cov_noise_1e-3_marginalization.png")
# savefig("rho1_posterior_cov_noise_1e-3_marginalization.pdf")

# savefig("rho1_posterior_cov_noise_1e-2_marginalization.png")
# savefig("rho1_posterior_cov_noise_1e-2_marginalization.pdf")

# savefig("rho1_posterior_cov_noise_1e-1_marginalization.png")
# savefig("rho1_posterior_cov_noise_1e-1_marginalization.pdf")

# savefig("rho1_posterior_cov_noise_1e0_marginalization.png")
# savefig("rho1_posterior_cov_noise_1e0_marginalization.pdf")

# savefig("rho1_posterior_cov_noise_1e1_marginalization.png")
# savefig("rho1_posterior_cov_noise_1e1_marginalization.pdf")

# savefig("rho1_posterior_cov_noise_1e2_marginalization.png")
# savefig("rho1_posterior_cov_noise_1e2_marginalization.pdf")


[norm(Γρ_I[:,:,i]) for i in 1:Nnoise]