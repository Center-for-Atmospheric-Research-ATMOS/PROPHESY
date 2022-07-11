## load the packages used in the estimation
# plotting
using PyPlot
rc("text", usetex=true)
rc("figure",max_open_warning=50)
using myPlot

# data manipulation (loading, writing, etc)
using Printf
using CSV
using DataFrames

# scientific package from the official Julia repositories
using LinearAlgebra
using StatsBase

# implemented scientific packages
using utilsFun  # for the softMax functions

# modeling XPS
using XPSpack

# tags
LOW_RES   = true               # set to true for computing the low resolution measurement models
MARG_UN   = (false & LOW_RES)    # set to true for computing the mean measurement operator as well as the covarainces
SHORT_RANGE = false              # select either wide range of attenuation lengths (false) or a restricted range more similar to experimental setup (true)
SIMULATE_DATA = true            # set true for simulating some data (several nose level)
SAVE_DATA = (true & SIMULATE_DATA)                # set to true for saving the generated variables
SAVE_MODEL = true               # set to true to save the models

MODEL_5   = false               # select the number of attenuation lengths probed
MODEL_10  = false
MODEL_20  = true

SAMPLE_MODEL = (true & LOW_RES) # if set to true, generate plenty of models by drawing attenuation lengths (only in low resolution, and does not compute the marginalization)
N_model_sample = 100;           # 100 should be more than enough for 5 attenuation length, but maybe not for 20

FLAG_0001 = (true & SIMULATE_DATA)               # selection of the profile (one must be true and the others false)
FLAG_0002 = (false & SIMULATE_DATA)
FLAG_0003 = (false & SIMULATE_DATA)
FLAG_0004 = (false & SIMULATE_DATA)

save_folder = "./";

# geometry setup
λe0 = 2.0e-3;        # reference penetration depth in μm
δr = 2.0e-3          # transition to vacuum layer thickness (let's set about 1 nm)
k0 = 10;             # compute the measurement model over a distance of k0*λe0
Nr = 201;            # number of discretization points in the radial dimension
Nr_lowres = 101 ;    # for the low resolution model
Nθ = 256;            # number of discretization points in the polar angle dimension
Ny = 256;            # number of discretization points in the cylinder axis dimension
μ0 = 10.0; # 20.0;   # radius of the cylinder
L = 100.0;           # height of the irradiated sample
x0 = sqrt(2.0)*100.0 # (x0,y0,z0) are the coordinates of the analyzer's apperture
y0 = 0.0;
z0 = 100.0;

# soon add the option of sphere or plane geometry?
save_folder = string(save_folder,"cylinder_radius_",μ0,"/")


# spacial discretization 
r = collect(range(μ0-k0*λe0,μ0+δr,length=Nr));
r_lowres = collect(range(μ0-k0*λe0,μ0+δr,length=Nr_lowres));
θ0 = atan(x0,z0)
θ = collect(range(θ0-π/2.0,θ0+π/2.0,Nθ));
y = collect(range(-L/2.0,L/2.0,length=Ny));

# concentration profiles (4 different cases)
ρ0 = 1.0
ρ_vac = 0.0
r_th  = 2.346; # 2.0


if FLAG_0001
    # ρA_1 = logistic.(1000.0reverse(μ0.-r).-r_th,ρ_vac,ρ0,2.0); # r_th=2.0
    ρA_1 = logistic.(1000.0reverse(δr.+μ0.-r).-r_th,ρ_vac,ρ0,2.0);
    exp_tag     = "0001"
end
if FLAG_0002
    # ρA_1 = logistic.(1000.0reverse(μ0.-r).-r_th,ρ_vac,ρ0,2.0) .+ 2.0exp.(-(1000.0reverse(μ0.-r).-1.0).^2. /(2.0*0.25^2)); # r_th=2.0
    ρA_1 = logistic.(1000.0reverse(δr.+μ0.-r).-r_th,ρ_vac,ρ0,2.0) .+ 2.0exp.(-(1000.0reverse(μ0.-r_surf).-0.0).^2. /(2.0*0.25^2));
    exp_tag     = "0002"
end
if FLAG_0003
    # ρA_1 = logistic.(1000.0reverse(μ0.-r).-r_th,ρ_vac,ρ0,2.0) .+ exp.(-(1000.0reverse(μ0.-r).-1.5).^2. /(2.0*0.5^2)); # r_th=2.0
    ρA_1 = logistic.(1000.0reverse(δr.+μ0.-r).-r_th,ρ_vac,ρ0,2.0) .+ exp.(-(1000.0reverse(μ0.-r_surf).-0.5).^2. /(2.0*0.5^2));
    exp_tag     = "0003"
end
if FLAG_0004
    # ρA_1 = exp.(-(1000.0reverse(μ0.-r).-2.5).^2. /(2.0*0.5^2));
    ρA_1 = exp.(-(1000.0(reverse(δr.+μ0.-r).-δr)).^2. /(2.0*0.5^2));
    exp_tag     = "0004"
end

# figure()
# plot(1000.0reverse(μ0.-r),ρA_1)

# measurement operator (only the geometrical term since the other comes as a multiplicative scalar estimated from the data)
if MODEL_5                                   # number of measurement (penetration depth)
    Ndata = 5;
end
if MODEL_10
    Ndata = 10;
end
if MODEL_20
    Ndata = 20;
end

if SHORT_RANGE                                          # bounds of the attenuation lengths
    λe1 = 1.3
    λe2 = 2.5
    save_folder = string(save_folder,"eal_",Ndata,"_restricted_range/")
else
    λe1 = 0.5;
    λe2 = 5.5; 
    save_folder = string(save_folder,"eal_",Ndata,"/")
end
λe = 1.0e-3collect(range(λe1,λe2,Ndata));              # attenuation length range




if LOW_RES
    δκ = [0.5; 1.0; 2.5]/100.0                             # relative error bound (attenuation length)
    H_lowres = zeros(Cdouble,Ndata,Nr_lowres);             # low resolution operator computed with corrupted attenuation length values
    save_folder_lowres = string(save_folder,"lowres/")     # base folder name where to save the models
    if SAVE_MODEL
        mkpath(save_folder_lowres)
        CSV.write(string(save_folder_lowres,"radial_discretization_lowres.csv"),DataFrame(reverse(r_lowres)',:auto);header=true)
    end
    if SAMPLE_MODEL
        for m in 1:min(N_model_sample,999)
            println(m,"/",min(N_model_sample,999)," model sample")
            for k in 1:length(δκ) # for each uncertainty level in the measurement model
                println(k,"/",length(δκ)," model uncertainty")
                for i in 1:Ndata
                    # draw the relative error 
                    κ = δκ[k]*(2.0*rand()-1.0);
                    # compute the operator with the corrupted attenuation length values
                    H_lowres[i,:],_,_,_,_ = cylinder_gain_H(r_lowres,θ,y,x0,y0,z0,μ0,(1.0+κ)*λe[i]);
                end
                global H_lowres = reverse(H_lowres,dims=2);
    
                # save the files where they belong
                if SAVE_MODEL
                    if (m<10)
                        mkpath(string(save_folder_lowres,"error_model_",δκ[k],"_percent/00",m,"/"))
                        CSV.write(string(save_folder_lowres,"error_model_",δκ[k],"_percent/00",m,"/","H_lowres.csv"),DataFrame(H_lowres,:auto);header=true)
                    elseif ((m>=10) & (m<100))
                        mkpath(string(save_folder_lowres,"error_model_",δκ[k],"_percent/0",m,"/"))
                        CSV.write(string(save_folder_lowres,"error_model_",δκ[k],"_percent/0",m,"/","H_lowres.csv"),DataFrame(H_lowres,:auto);header=true)
                    else
                        mkpath(string(save_folder_lowres,"error_model_",δκ[k],"_percent/",m,"/"))
                        CSV.write(string(save_folder_lowres,"error_model_",δκ[k],"_percent/",m,"/","H_lowres.csv"),DataFrame(H_lowres,:auto);header=true)
                    end
                end
            end
        end
    else
        Nλ = 21;                                               # number of discretization point for the distribution of each attenuation length
        μH = Array{Cdouble}(undef,Nr_lowres,Ndata);            # mean value operator (w.r.t. λe distribution)
        ΓH = Array{Cdouble}(undef,Nr_lowres,Nr_lowres,Ndata);  # covariance of the corrupted measurement operator (w.r.t. λe distribution)
        for k in 1:length(δκ) # for each uncertainty level in the measurement model
            for i in 1:Ndata
                # draw the relative error 
                κ = δκ[k]*(2.0*rand()-1.0);
                # compute the operator with the corrupted attenuation length values
                H_lowres[i,:],_,_,_,_ = cylinder_gain_H(r_lowres,θ,y,x0,y0,z0,μ0,(1.0+κ)*λe[i]);
                if MARG_UN
                    # compute the mean operator and covariance for the marginalization of error (to make it realistic, I use the corrupted values as the center of the distribution instead of the true values which would not be accessible)
                    local λrange = (1.0+κ)*collect(range((1.0-δκ[k])*λe[i],(1.0+δκ[k])*λe[i],length=Nλ))
                    local Pλ = (1.0/((λrange[2]-λrange[1])*Nλ))*ones(Cdouble,Nλ)
                    ΓH[:,:,i], μH[:,i]   = cov_H_cylinder(r_lowres, θ,y,x0,y0,z0,μ0,λrange,Pλ)
                end
            end
            global H_lowres = reverse(H_lowres,dims=2);
            if MARG_UN
                global ΓH = reverse(ΓH,dims=(1,2));
                global μH = reverse(μH,dims=1);
            end

            # save the files where they belong
            if SAVE_MODEL
                mkpath(string(save_folder_lowres,"error_model_",δκ[k],"_percent/"))
                CSV.write(string(save_folder_lowres,"error_model_",δκ[k],"_percent/","H_lowres.csv"),DataFrame(H_lowres,:auto);header=true)
                if MARG_UN
                    mkpath(string(save_folder_lowres,"error_model_",δκ[k],"_percent/cov/"))
                    CSV.write(string(save_folder_lowres,"error_model_",δκ[k],"_percent/","mean_H_lowres.csv"),DataFrame(μH,:auto);header=true)
                    for i in 1:Ndata
                        CSV.write(string(save_folder_lowres,"error_model_",δκ[k],"_percent/cov/","cov_H_lowres_",i,".csv"),DataFrame(ΓH[:,:,i],:auto);header=false)
                    end
                end
            end
        end
    end
end





if SIMULATE_DATA
    H_highres = zeros(Cdouble,Ndata,Nr);                   # high resolution operator used for the data simulation (the attenuation length used for the computation of this operator are the reference ones)
    for i in 1:Ndata
        # high resolution measurement operator (using true λe values)
        H_highres[i,:],_,_,_,_ = cylinder_gain_H(r,θ,y,x0,y0,z0,μ0,λe[i]);
    end
    H_highres = reverse(H_highres,dims=2);
    if SAVE_MODEL
        mkpath(save_folder)
        CSV.write(string(save_folder,"radial_discretization.csv"),DataFrame(reverse(r)',:auto);header=true)
        CSV.write(string(save_folder,"attenuation_length.csv"),DataFrame(λe',:auto);header=false)
        CSV.write(string(save_folder,"H_highres.csv"),DataFrame(H_highres,:auto);header=true)
    end
    save_folder_data = string(save_folder,exp_tag,"/")
    # generate some data (data point and covariance)
    Nnoise = 200;
    SNRmodulation = 1.0 .+ (λe.-λe[1])./(λe[end]-λe[1]);
    SNR_level = [5.0; 10.0; 50.0; 100.0; 500.0; 1000.0];
    Signal_level = abs.(H_highres*ρA_1);

    
    for k in 1:length(SNR_level)
        # figure()
        y_data = zeros(Cdouble,Nnoise,Ndata);
        ΓIsqrt = diagm(Signal_level./(SNR_level[k]*SNRmodulation));
        y_data[1,:] = H_highres*ρA_1;
        for i in 2:Nnoise
            y_data[i,:] = y_data[1,:] + ΓIsqrt*randn(Cdouble,Ndata);
            # scatter(collect(1:Ndata),y_data[i,:])
        end
        y_data[y_data.<0.0] = -y_data[y_data.<0.0]; # just making sure the integrated peak area are positive values
        if SAVE_DATA
            mkpath(string(save_folder_data,"/noise_level_",σ_level[k],"/"))
            CSV.write(string(save_folder_data,"/noise_level_",σ_level[k],"/","repeated_data.csv"),DataFrame([λe';y_data],:auto);header=true)
        end
    end
    if SAVE_DATA
        CSV.write(string(save_folder_data,"concentration_profile.csv"),DataFrame(ρA_1',:auto);header=true)
    end
end


# 
# spectrum generation
# 

dKe = 0.05;
Be = collect(286.0:dKe:298.0);
μBe = [290.2; 292.0; 293.5] # assume that the peaks are at the same location for every photon energy

σ_be = [0.65; 0.45; 0.6];

p_peak = zeros(Cdouble,Ndata,3);
p_peak[:,1] = 0.62 .+ (0.4-0.62)*(λe.-λe[1])./(λe[end]-λe[1]);
p_peak[:,2] = 0.35 .+ (0.3-0.35)*(λe.-λe[1])./(λe[end]-λe[1]);
p_peak[:,3] = 1.0 .- (p_peak[:,1]+p_peak[:,2]);

figure()
for i in 1:Ndata
    σ_peak_1 = (1.0/sqrt(2.0π*σ_be[1]^2))*exp.(-0.5*(Be.-μBe[1]).^2/(2.0σ_be[1]^2));
    σ_peak_2 = (1.0/sqrt(2.0π*σ_be[2]^2))*exp.(-0.5*(Be.-μBe[2]).^2/(2.0σ_be[2]^2));
    σ_peak_3 = (1.0/sqrt(2.0π*σ_be[3]^2))*exp.(-0.5*(Be.-μBe[3]).^2/(2.0σ_be[3]^2));
    plot(Be,p_peak[i,1]*σ_peak_1+p_peak[i,2]*σ_peak_2+p_peak[i,3]*σ_peak_3)
end

