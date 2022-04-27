## load the packages used in the estimation
# plotting
using PyPlot
fm = PyPlot.matplotlib.font_manager.json_load("/home/matthew/.cache/matplotlib/fontlist-v310.json") # TODO: look for the path automatically
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
using XPSinv


# geometry setup
λe0 = 2.0e-3;        # reference penetration depth in μm
k0 = 10;             # compute the measurement model over a distance of k0*λe0
Nr = 201;            # number of discretization points in the radial dimension
Nr_lowres = 101 ;    # for the low resolution model
Nθ = 256;            # number of discretization points in the polar angle dimension
Ny = 256;            # number of discretization points in the cylinder axis dimension
μ0 = 20.0            # radius of the cylinder
L = 200.0;           # height of the irradiated sample
x0 = sqrt(2.0)*100.0 # (x0,y0,z0) are the coordinates of the analyzer's apperture
y0 = 0.0;
z0 = 100.0

# spacial discretization 
r = collect(range(μ0-k0*λe0,μ0,length=Nr))
r_lowres = collect(range(μ0-k0*λe0,μ0,length=Nr_lowres));
θ0 = atan(x0,z0)
θ = collect(range(θ0-π/2.0,θ0+π/2.0,Nθ));
y = collect(range(-L/2.0,L/2.0,length=Ny));

# concentration profiles (4 different cases)
ρA_1 = logistic.(1000.0reverse(μ0.-r).-2.0,0.0,1.0,2.0);
# ρA_1 = logistic.(1000.0reverse(μ0.-r).-2.0,0.0,1.0,2.0) .+ 2.0exp.(-(1000.0reverse(μ0.-r).-1.0).^2. /(2.0*0.25^2));
# ρA_1 = logistic.(1000.0reverse(μ0.-r).-2.0,0.0,1.0,2.0) .+ exp.(-(1000.0reverse(μ0.-r).-1.5).^2. /(2.0*0.5^2));
# ρA_1 = exp.(-(1000.0reverse(μ0.-r).-2.5).^2. /(2.0*0.5^2));


# measurement operator (only the geometrical term since the other comes as a multiplicative scalar estimated from the data)
Ndata = 20 # 5 # 10                                    # number of measurement (penetration depth)
λe1 = 0.5; # 1.3                                       # bounds of the attenuation lengths
λe2 = 5.5; # 2.5
λe = 1.0e-3collect(range(λe1,λe2,Ndata));              # attenuation length range
Nλ = 21;                                               # number of discretization point for the distribution of each attenuation length
δκ = 2.5/100.0; # 1.0/100.0; # 0.5/100.0; #            # relative error bound (attenuation length)

H_highres = zeros(Cdouble,Ndata,Nr);                   # high resolution operator used for the data simulation (the attenuation length used for the computation of this operator are the reference ones)
H_lowres = zeros(Cdouble,Ndata,Nr_lowres);             # low resolution operator computed with corrupted attenuation length values
μH = Array{Cdouble}(undef,Nr_lowres,Ndata);            # mean value operator (w.r.t. λe distribution)
ΓH = Array{Cdouble}(undef,Nr_lowres,Nr_lowres,Ndata);  # covariance of the corrupted measurement operator (w.r.t. λe distribution)
for i in 1:Ndata
    # high resolution measurement operator (using true λe values)
    H_highres[i,:],_,_,_,_ = cylinder_gain_H(r,θ,y,x0,y0,z0,μ0,λe[i]);
    # draw the relative error 
    κ = δκ*(2.0*rand()-1.0);
    # compute the operator with the corrupted attenuation length values
    H_lowres[i,:],_,_,_,_ = cylinder_gain_H(r_lowres,θ,y,x0,y0,z0,μ0,(1.0+κ)*λe[i]);
    # compute the mean operator and covariance for the marginalization of error (to make it realistic, I use the corrupted values as the center of the distribution instead of the true values which would not be accessible)
    local λrange = (1.0+κ)*collect(range((1.0-δκ)*λe[i],(1.0+δκ)*λe[i],length=Nλ))
    local Pλ = (1.0/((λrange[2]-λrange[1])*Nλ))*ones(Cdouble,Nλ)
    ΓH[:,:,i], μH[:,i]   = cov_H_cylinder(r_lowres, θ,y,x0,y0,z0,μ0,λrange,Pλ)
end
H_highres = reverse(H_highres,dims=2);
H_lowres = reverse(H_lowres,dims=2);
ΓH = reverse(ΓH,dims=(1,2));
μH = reverse(μH,dims=1);


# CSV.write("radial_discretization.csv",DataFrame(reverse(r)',:auto);header=true)
# CSV.write("radial_discretization_lowres.csv",DataFrame(reverse(r_lowres)',:auto);header=true)
# CSV.write("attenuation_length.csv",DataFrame(λe',:auto);header=false)
# CSV.write("H_highres.csv",DataFrame(H_highres,:auto);header=true)

# CSV.write("H_lowres_25.csv",DataFrame(H_lowres,:auto);header=true)
# CSV.write("mean_H_lowres_25.csv",DataFrame(μH,:auto);header=true)
# for i in 1:Ndata
#     CSV.write(string("cov_H_lowres_25_",i,".csv"),DataFrame(ΓH[:,:,i],:auto);header=false)
# end

# CSV.write("H_lowres_1.csv",DataFrame(H_lowres,:auto);header=true)
# CSV.write("mean_H_lowres_1.csv",DataFrame(μH,:auto);header=true)
# for i in 1:Ndata
#     CSV.write(string("cov_H_lowres_1_",i,".csv"),DataFrame(ΓH[:,:,i],:auto);header=false)
# end

# CSV.write("H_lowres_05.csv",DataFrame(H_lowres,:auto);header=true)
# CSV.write("mean_H_lowres_05.csv",DataFrame(μH,:auto);header=true)
# for i in 1:Ndata
#     CSV.write(string("cov_H_lowres_05_",i,".csv"),DataFrame(ΓH[:,:,i],:auto);header=false)
# end


# generate some data (data point and covariance)
Nnoise = 200;
# σnoise = 0.001*ones(Cdouble,Nnoise);
σnoise = 0.005*ones(Cdouble,Nnoise);
# σnoise = 0.01*ones(Cdouble,Nnoise);
# σnoise = 0.05*ones(Cdouble,Nnoise);
# σnoise = 0.1*ones(Cdouble,Nnoise);

y_data = zeros(Cdouble,Nnoise,Ndata);
ΓI = zeros(Cdouble,Ndata,Ndata,Nnoise);
ΓIsqrt = zeros(Cdouble,Ndata,Ndata,Nnoise);
# ΓIinv = zeros(Cdouble,Ndata,Ndata,Nnoise);
for i in 1:Nnoise
    ΓI[:,:,i] = σnoise[i]^2*diagm(ones(Cdouble,Ndata)); # iid data noise
    ΓIsqrt[:,:,i] = sqrt(ΓI[:,:,i]);
    # ΓIinv[:,:,i] = inv(ΓI[:,:,i]);
    if i==1
        y_data[i,:] = H_highres*ρA_1;
    else
        y_data[i,:] = H_highres*ρA_1 + ΓIsqrt[:,:,i]*randn(Cdouble,Ndata);
    end
end
y_data[y_data.<0.0] = -y_data[y_data.<0.0];


# CSV.write("repeated_data.csv",DataFrame(y_data,:auto);header=true)
# CSV.write("noise_level.csv",DataFrame(σnoise',:auto);header=false)
