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
# using StatsBase

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


LH = cholesky(ΓH_pos)

sqrtΓH = sqrt(ΓH);

figure();
plot(r,μH.+real(sqrtΓH)*randn(Cdouble,Nr,20))


# norm(real(sqrtΓH))
# norm(imag(sqrtΓH))

# ℓ = 1;
# Rℓλ = exp.(-D[ℓ,:,:])