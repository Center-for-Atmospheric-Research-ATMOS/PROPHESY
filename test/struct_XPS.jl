## load the packages used in the estimation
# plotting
using PyPlot
fm = PyPlot.matplotlib.font_manager.json_load("/home/matthew/.cache/matplotlib/fontlist-v310.json") # TODO: look for the path automatically
# fm = PyPlot.matplotlib.font_manager.json_load("/home/mattoz/.cache/matplotlib/fontlist-v310.json")
# fm.findfont("serif", rebuild_if_missing=false)
# fm.findfont("serif", fontext="afm", rebuild_if_missing=false)
rc("font",family="serif",serif="Computer Modern Roman")
rc("text", usetex=true)
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

# Dict{Int64,XPSsetup}();


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
x0 = sqrt(2.0)*100.0#μ0
y0 = 0.0;
z0 = 100.0#μ0
r = collect(range(μ0-k0*λe0,μ0,length=Nr))
θ0 = atan(x0,z0)
θ = collect(range(θ0-π/2.0,θ0+π/2.0,Nθ));
y = collect(range(-L/2.0,L/2.0,length=Ny));

wsGeom = cylinderGeom(x0,y0,z0,μ0,r,θ,y)

# TODO: create a model from the given elements
Hr,Hrθy,Arn,Aθj,Ayk = cylinder_gain_H(r,θ,y,x0,y0,z0,μ0,λe0);



fig,ax,pcm,cax,cb = imshowDataPolar(1,r,θ,Hrθy[:,:,128];cb_ax_loc=(0.25, .37));
# ax.set_rticks([99.97, 99.98, 99.99, 100.0])
ax.set_ylim(99.97,100.0)



function acqModel(wsAcq::XPSacq,wsGeom::cylinderGeom)
    Hr,Hrθy,_,_,_ = cylinder_gain_H(wsGeom.r,wsGeom.θ,wsGeom.y,wsGeom.x0,wsGeom.y0,wsGeom.z0,wsGeom.μ0,wsAcq.λe);
    wsAcq.T*wsAcq.α*wsAcq.Fν*σν*Hr'
end


# a few things about model uncertainty

# """
#     cov_H_cylinder()

#     computes the covariance matrix of the geometrical structure...
#     A = ∭ ρ(r,θ,y) e^{\frac{d_P(r,θ,y)}{λ}} r dr dθ dy ≃ Hρ
#     where H = [H_1 H_2 … H_Nr]† and 
#     H_n(λ) = r_n ∑_j ∑_k e^{\frac{d_P(r_n,θ_j,y_k)}{λ}} ∭ e_n(r) e_j(θ) e_k(y) dr dθ dy
#     ΓH = cov(H) = \mathbb{E} [H×H†] - \mathbb{E}[H]×\mathbb{E} [H†]
# """
function cov_H_cylinder(r::Array{Cdouble,1},θ::Array{Cdouble,1},y::Array{Cdouble,1},x0::Cdouble,y0::Cdouble,z0::Cdouble,μ0::Cdouble,λe0::Cdouble)
    # the distance for each point of the space discretization
    Nr = length(r);
    Nθ = length(θ);
    Ny = length(y);
    D = zeros(Cdouble,Nr,Nθ,Ny);
    for k in 1:Ny
        D[:,:,k] = d_cylinder_P(r,θ,y[k],x0,y0,z0,μ0);
    end
    Arn = 0.5*r.*[r[2]-r[1]; r[3:end]-r[1:end-2]; r[end]-r[end-1]]; # rdr
    Aθj = 0.5*[θ[2]-θ[1]; θ[3:end]-θ[1:end-2]; θ[end]-θ[end-1]];    # dθ
    Ayk = 0.5*[y[2]-y[1]; y[3:end]-y[1:end-2]; y[end]-y[end-1]];    # dy

    # attenuation length distribution
    Nλ = 21;
    λ = collect(range(0.9λe0,1.1λe0,length=Nλ)); # should be given as an argument because 
    Aλ = 0.5*[λ[2]-λ[1]; λ[3:end]-λ[1:end-2]; λ[end]-λ[end-1]];
    Pλ = (1.0/sum(Aλ))ones(Cdouble,Nλ);          # should be given as an argument

    # compute the integration operator for the discretized attenuation space
    H = zeros(Cdouble,Nr,Nλ);
    Djk = Aθj*Ayk'; # integration over θ and y 
    for m in 1:Nr
        for s in 1:Nλ
            H[m,s] = Arn[m]*sum(Djk.*exp.(-D[m,:,:]/λ[s]))
        end
    end

    # mean operator 
    μH = H*(Pλ.*Aλ);

    # square
    HHt = zeros(Cdouble,Nr,Nr);
    for l in 1:Nr
        HHt[l,l] = (H[l,:].^2 .*Pλ)'*Aλ
        for m in l+1:Nr
            HHt[l,m] = (H[l,:].*H[m,:].*Pλ)'*Aλ
            HHt[m,l] = HHt[l,m]
        end
    end

    # return the covariance
    HHt - μH*μH', μH
end



ΓH,μH = cov_H_cylinder(r,θ,y,x0,y0,z0,μ0,λe0);

figure()
imshow(ΓH)
colorbar()

evals = eigvals(ΓH)

figure()
semilogy(abs.(evals))

svals = svdvals(ΓH)

figure()
semilogy(svals)

figure();
semilogy(r,μH,color="tab:blue")
fill_between(r,μH-sqrt.(diag(ΓH)),μH+sqrt.(diag(ΓH)),alpha=0.5,color="tab:blue",label="uncertainty")

FΓH = eigen(ΓH);

pos_val = zeros(Cdouble,Nr); # = FΓH.values[]
th_val = 1.0e-5*maximum(FΓH.values);
pos_val[FΓH.values.>th_val] .= FΓH.values[FΓH.values.>th_val];
pos_val[FΓH.values.<=th_val] .= th_val;

pos_val[FΓH.values.<=0.0] .= 1.0e-20; # minimum(FΓH.values[FΓH.values.>0])*0.5*rand(Cdouble,Nr-length(FΓH.values[FΓH.values.>0]));

ΓH_pos = FΓH.vectors*diagm(pos_val)*FΓH.vectors';
# ΓH_pos = 0.5*(ΓH_pos+ΓH_pos');

figure(); imshow(ΓH); colorbar()
figure(); imshow(ΓH_pos); colorbar()
figure(); imshow(abs.(ΓH-ΓH_pos)); colorbar()


# LH = cholesky(ΓH_pos)

sqrtΓH = sqrt(ΓH);

figure();
plot(r,μH.+real(sqrtΓH)*randn(Cdouble,Nr,20))


# norm(real(sqrtΓH))
# norm(imag(sqrtΓH))

# ℓ = 1;
# Rℓλ = exp.(-D[ℓ,:,:])


##
## OK, going in a different direction for estimating the uncertainty (posterior covariance): sampling the a priori with rejection mechanism based on likelihood (or a posteriori)
##


# simulate some data (one point in the kinetic energy spectrum for four different concentration profiles)
ρA_1 = logistic.(1000.0reverse(μ0.-r).-2.0,0.0,1.0,2.0);
ρA_2 = logistic.(1000.0reverse(μ0.-r).-2.0,0.0,1.0,2.0) .+ 2.0exp.(-(1000.0reverse(μ0.-r).-1.0).^2. /(2.0*0.25^2));
ρA_3 = logistic.(1000.0reverse(μ0.-r).-2.0,0.0,1.0,2.0) .+ exp.(-(1000.0reverse(μ0.-r).-1.5).^2. /(2.0*0.5^2));
ρA_4 = exp.(-(1000.0reverse(μ0.-r).-2.5).^2. /(2.0*0.5^2));

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


# Lprior = cholesky(Γprior)

Dprior = D2nd(Nr) # D2nd(Nr+2)[:,2:end-1];
# Bprior = 1.0e8Dprior'*inv(Γprior[2:end-1,2:end-1])*Dprior; # 1.0e-8
Bprior = 1.0e-16Dprior'*inv(Γprior[2:end-1,2:end-1])*Dprior; # 1.0e-8
Cprior = inv(Bprior);
Dsqrt = real(sqrt(Cprior));

Ndata = 5 # 25
H_dummy = zeros(Cdouble,Ndata,Nr);
H_better = zeros(Cdouble,Ndata,Nr);
Γbetter  = zeros(Cdouble,Ndata,Nr,Nr);
λbetter  = 1.0e-3*[1.0; 1.5; 2.0; 2.5; 3.0];
for i in 1:Ndata
    # H_dummy[i,:] = 5.0exp.(-collect(range(0.0,Nr,length=Nr))./(0.5*i))
    # H_dummy[i,:] = 5.0exp.(-0.5*(collect(range(0.0,Nr,length=Nr)).-0.5i.-15.0).^2 ./(5.0^2))
    # H_dummy[i,:] = 1.0(exp.(-collect(range(0.0,Nr,length=Nr))./(0.5*(i+1)+1.0))-exp.(-collect(range(0.0,Nr,length=Nr))./(0.5*i+1.0)))
    H_dummy[i,:] = exp.(-collect(range(0.0,Nr,length=Nr))./(0.5*i+2.0))
    
    Γbetter[i,:,:],H_better[i,:] = cov_H_cylinder(r,θ,y,x0,y0,z0,μ0,λbetter[i]);
end


figure(); imshow(H_dummy); colorbar()

figure(); plot(μ0.-r,H_dummy')

H_better = reverse(H_better,dims=2);
figure(); plot(μ0.-r,reverse(H_better,dims=2)')

Γbetter = reverse(Γbetter,dims=(2,3));

for i in 1:Ndata
    figure(); imshow(Γbetter[i,:,:]); colorbar()
end


ΓI = (2.0e-2^2)*diagm(ones(Cdouble,Ndata)); # iid data noise
# ΓI = (1.1^2)*diagm(ones(Cdouble,Ndata));
ΓIsqrt = sqrt(ΓI);
detΓI = det(ΓI);
ΓIinv = inv(ΓI);

if false
    y_data_1 = H_dummy*ρA_1 + ΓIsqrt*randn(Cdouble,Ndata);
    y_data_1[y_data_1.<0.0] = -y_data_1[y_data_1.<0.0];
else
    y_data_1 = H_better*ρA_1 + ΓIsqrt*randn(Cdouble,Ndata);
    y_data_1[y_data_1.<0.0] = -y_data_1[y_data_1.<0.0];
end


figure(); plot(y_data_1)

# ρ_all = zeros(Cdouble,Ns+1,Nr);
# square root matrix of the generative covariance matrix (the covariance in the distribution used for generating new samples)
σw = 10.e-2 # 0.001; # small compared with the amplitude of the state 
w = σw*ones(Cdouble,Nr); # not optimal because we know that the concentration varies more in the region near the surface rather than deep in the sample
# w = σw*(1.0.-ρA_1.+0.1); # too optimal because we known the solution (but the general sigmoid shape could be used because it's not a big a priori)
Γsqrt = real(sqrt(corrCovariance(w;cor_len=15.0)));
p0 = 0.5 # 0.02; #starting acceptance rate of uphill moves
ρB = [ρA_1[1]; ρA_1[end]];
σB = [0.01; 0.01];
Ns = 1000000;
if false
    ρ_all = samplePosterior(0.0ρA_1,Γsqrt,p0*ones(Cdouble,Ns),y_data_1,ΓIinv,H_dummy,Bprior,ρB,σB;Ns=Ns,psmooth=1.999);
else
    ρ_all = samplePosteriorModelMargin(0.0ρA_1,Γsqrt,p0*ones(Cdouble,Ns),y_data_1,ΓIinv,H_better,Γbetter,Bprior,ρB,σB;Ns=Ns,psmooth=1.999);
end
# ρ_all = samplePosterior(ρA_1,Γsqrt,p0,y_data_1,ΓIinv,H_dummy,Bprior,ρB,σB;Ns=Ns,psmooth=1.999);

μρ_IG = dropdims(mean(ρ_all,dims=1),dims=1);
Γρ_IG = cov(ρ_all.-μρ_IG');

# figure(); plot(μ0.-r,ρ_all[end-20:end,:]')

figure();
plot(μ0.-r,ρA_1)
plot(μ0.-r,μρ_IG,color="blue")
fill_between(μ0.-r,ρA_1-sqrt.(diag(Γρ_IG)),ρA_1+sqrt.(diag(Γρ_IG)),alpha=0.5,color="tab:blue",label="uncertainty")
# fill_between(μ0.-r,μρ_IG-sqrt.(diag(Γρ_IG)),μρ_IG+sqrt.(diag(Γρ_IG)),alpha=0.5,color="tab:blue",label="uncertainty")

figure(); imshow(Γρ_IG); colorbar()

figure(); plot(μ0.-r,sqrt.(diag(Γρ_IG)))


# TODO: observe the burn in period and don't use it for the computation of the mean and the covariance
Elikelihood  = zeros(Cdouble,Ns+1);
EpriorSmooth = zeros(Cdouble,Ns+1);
EpriorVal    = zeros(Cdouble,Ns+1);
EpriorModel  = zeros(Cdouble,Ns+1);
for i in 1:Ns+1
    Elikelihood[i]  = (y_data_1-H_dummy*ρ_all[i,:])'*ΓIinv*(y_data_1-H_dummy*ρ_all[i,:])
    EpriorSmooth[i] = ρ_all[i,:]'Bprior*ρ_all[i,:]
    EpriorVal[i] = ((ρ_all[i,1]  -ρB[1])^2)/(σB[1]^2) + ((ρ_all[i,end]  -ρB[2])^2)/(σB[2]^2)
    for g in 1:Ndata
        EpriorModel[i] = EpriorModel[i] + ΓIinv[i,i]*(ρ_all[i,:]'*Γbetter[i,:,:]*ρ_all[i,:])
    end
end

Etot = Elikelihood+EpriorSmooth+EpriorVal;
maxProbDens = maximum(exp.(-0.5*(Elikelihood+EpriorSmooth+EpriorVal)))

idx_best = findall(exp.(-0.5*(Elikelihood+EpriorSmooth+EpriorVal)).>=0.1maxProbDens)

figure()
semilogx(collect(1:Ns+1),Elikelihood,label="likelihood")
semilogx(collect(1:Ns+1),EpriorSmooth,label="smoothness a priori")
semilogx(collect(1:Ns+1),EpriorVal,label="values a priori")
semilogx(collect(1:Ns+1),EpriorModel,label="values a priori")
semilogx(collect(1:Ns+1),Elikelihood+EpriorSmooth+EpriorVal+EpriorModel,label="total")
legend()

# figure()
# plot(Elikelihood[idx_best],label="likelihood")
# plot(EpriorSmooth[idx_best],label="smoothness a priori")
# plot(EpriorVal[idx_best],label="values a priori")
# plot(Elikelihood[idx_best]+EpriorSmooth[idx_best]+EpriorVal[idx_best],label="total")
# legend()

# try and plot the probability: not a good idea because the values in the exp are often to small, returning a lot of 0s
figure()
plot(exp.(-0.5*(Elikelihood+EpriorSmooth+EpriorVal))) 


mean(Elikelihood[100000:end])
mean(EpriorSmooth[100000:end])
mean(EpriorVal[100000:end])

μρ_IG = dropdims(mean(ρ_all[100000:end,:],dims=1),dims=1);
Γρ_IG = cov(ρ_all[100000:end,:]);


figure();
plot(μ0.-r,ρA_1)
plot(μ0.-r,μρ_IG,color="blue")
fill_between(μ0.-r,ρA_1-sqrt.(diag(Γρ_IG)),ρA_1+sqrt.(diag(Γρ_IG)),alpha=0.5,color="tab:blue",label="uncertainty")
# fill_between(μ0.-r,μρ_IG-sqrt.(diag(Γρ_IG)),μρ_IG+sqrt.(diag(Γρ_IG)),alpha=0.5,color="tab:blue",label="uncertainty")

figure(); imshow(Γρ_IG); colorbar()

figure(); plot(μ0.-r,sqrt.(diag(Γρ_IG)))


# acceptance rate over the uphill moves

dE = Etot[2:end]-Etot[1:end-1];
τ = (length(dE[dE.>=0.0])-length(dE[dE.==0.0]))/(length(dE[dE.>=0.0]))

# overall acceptance rate over the iterations
Eb = 1.0*(dE.>0.0)

Eτ = zeros(Cdouble,Ns);
for i in 1:Ns-999
    Eτ[i] = mean(Eb[i:i+999])
end

figure(); plot(collect(1.0:Ns),Eτ) # looks quite good, the burning in period seems really negligible in this case

