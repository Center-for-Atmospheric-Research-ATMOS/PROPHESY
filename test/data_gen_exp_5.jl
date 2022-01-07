## load the packages used in the estimation
# plotting
using PyPlot
# set the font in the plots
# PyPlot.matplotlib.rcParams["font.family"]="serif"
# rc("font", family="serif")
# PyPlot.matplotlib.rcParams["mathtext.fontset"]="Computer Modern Roman"
# rc("font", serif="Computer Modern Roman")
# PyPlot.matplotlib.rcParams["text.usetex"]
# fm = matplotlib.font_manager.json_load("/path/to/fontlist-v300.json")
# fm.findfont("serif", rebuild_if_missing=False)
# fm.findfont("serif", fontext="afm", rebuild_if_missing=False)
fm = PyPlot.matplotlib.font_manager.json_load("/home/mattoz/.cache/matplotlib/fontlist-v310.json")
fm.findfont("serif", rebuild_if_missing=false)
fm.findfont("serif", fontext="afm", rebuild_if_missing=false)
rc("font",family="serif",serif="Computer Modern Roman")
rc("text", usetex=true)
using myPlot

# data manipulation (loading, writing, etc)
using CSV
using Dates
using DataFrames # for loading and saving the results, e.g. CSV.write("measured_psd.csv",DataFrame(PSD_meas'); writeheader=false)
using Printf

# scientific package from the official Julia repositories
using LinearAlgebra
using Statistics
using DSP
using SpecialMatrices
using Polynomials
using StatsBase

# implemented scientific packages
using utilsFun  # for the softMax functions

# modeling XPS
using XPSpack

# what to save
SAVE_MODEL = false;
SOME_INV   = false; #WARNING: not the peak area version, requires XPSinv (not yet committed)
PLOT_FIG   = false


d_c1s_Ek = CSV.File(string("data/","C1s on Ek, Na-prop and Prop.csv"); delim=",", header=true) |> DataFrame
Ny = length(d_c1s_Ek[:,1]);
D_2nd = diagm(Ny-2,Ny,1 => 2ones(Cdouble,Ny-2), 0 => -ones(Cdouble,Ny-2) ,2 => -ones(Cdouble,Ny-2)); # D2nd(Ny);
κb = 0.01;
λb = 1.0e5;
z_baseline_1 = baseline_removal(d_c1s_Ek[:,2],λb,κb,D_2nd);
if PLOT_FIG
   figure();
   scatter(d_c1s_Ek[:,1], d_c1s_Ek[:,2])
   plot(d_c1s_Ek[:,1], d_c1s_Ek[:,2]-z_baseline_1)
   plot(d_c1s_Ek[:,1], z_baseline_1)
end

# κb = 0.01;
# λb = 1.0e5; # 1.0e3
z_baseline_2 = baseline_removal(d_c1s_Ek[:,4],λb,κb,D_2nd);
if PLOT_FIG
   figure();
   scatter(d_c1s_Ek[:,3], d_c1s_Ek[:,4])
   plot(d_c1s_Ek[:,3], d_c1s_Ek[:,4]-z_baseline_2)
   plot(d_c1s_Ek[:,3], z_baseline_2)
end

# κb = 0.01;
# λb = 1.0e5;
z_baseline_3 = baseline_removal(d_c1s_Ek[:,6],λb,κb,D_2nd);
if PLOT_FIG
   figure();
   scatter(d_c1s_Ek[:,5], d_c1s_Ek[:,6])
   plot(d_c1s_Ek[:,5], d_c1s_Ek[:,6]-z_baseline_3)
   plot(d_c1s_Ek[:,5], z_baseline_3)
end

# κb = 0.01;
# λb = 1.0e5; # 1.0e4
z_baseline_4 = baseline_removal(d_c1s_Ek[:,8],λb,κb,D_2nd);
if PLOT_FIG
   figure();
   scatter(d_c1s_Ek[:,7], d_c1s_Ek[:,8])
   plot(d_c1s_Ek[:,7], d_c1s_Ek[:,8]-z_baseline_4)
   plot(d_c1s_Ek[:,7], z_baseline_4)
end

# κb = 0.01;
# λb = 1.0e5;
z_baseline_5 = baseline_removal(d_c1s_Ek[:,10],λb,κb,D_2nd);
if PLOT_FIG
   figure();
   scatter(d_c1s_Ek[:,9], d_c1s_Ek[:,10])
   plot(d_c1s_Ek[:,9], d_c1s_Ek[:,10]-z_baseline_5)
   plot(d_c1s_Ek[:,9], z_baseline_5)
end

x_all = [reverse(d_c1s_Ek[:,1]).+60.0; reverse(d_c1s_Ek[:,3]).+210.0; reverse(d_c1s_Ek[:,5]).+410.0; reverse(d_c1s_Ek[:,7]).+610.0; reverse(d_c1s_Ek[:,9]).+910.0];
y_all_with_baseline = [reverse(d_c1s_Ek[:,2]); reverse(d_c1s_Ek[:,4]); reverse(d_c1s_Ek[:,6]); reverse(d_c1s_Ek[:,8]); reverse(d_c1s_Ek[:,10])];
z_baseline = [reverse(z_baseline_1); reverse(z_baseline_2); reverse(z_baseline_3); reverse(z_baseline_4); reverse(z_baseline_5)];
y_all = y_all_with_baseline - z_baseline;
if PLOT_FIG
   figure();
   scatter(x_all, y_all_with_baseline)
   plot(x_all,z_baseline)
   plot(x_all, y_all_with_baseline-z_baseline)
   plot(x_all, z_baseline)
   plot([x_all[1]; x_all[end]],zeros(Cdouble,2))
end

τm = [1.0/3.0; 1.0/3.0; 1.0/3.0];
μm = [290.3; 291.9; 293.5];
σm = [290.3; 291.9; 293.5]/500.0;

NN = length(d_c1s_Ek[:,1]);
τt = zeros(Cdouble,5,3)
μt = zeros(Cdouble,5,3)
σt = zeros(Cdouble,5,3)
for i in 1:5
   # estimate the peaks centers and spreads
   τt[i,:],μt[i,:],σt[i,:] = EM_peaks(d_c1s_Ek[:,(2i)-1],d_c1s_Ek[:,2i]-reverse(z_baseline[(i-1)*NN+1:i*NN]),τm,μm,σm,100)
   # plt
   x_dist_1 = (τt[i,1]/σt[i,1])*exp.(-0.5*((d_c1s_Ek[:,(2i)-1].-μt[i,1])/σt[i,1]).^2);
   x_dist_2 = (τt[i,2]/σt[i,2])*exp.(-0.5*((d_c1s_Ek[:,(2i)-1].-μt[i,2])/σt[i,2]).^2);
   x_dist_3 = (τt[i,3]/σt[i,3])*exp.(-0.5*((d_c1s_Ek[:,(2i)-1].-μt[i,3])/σt[i,3]).^2);
   if PLOT_FIG
      figure();
      plot(d_c1s_Ek[:,(2i)-1],x_dist_1)
      plot(d_c1s_Ek[:,(2i)-1],x_dist_2)
      plot(d_c1s_Ek[:,(2i)-1],x_dist_3)
      plot(d_c1s_Ek[:,(2i)-1],x_dist_1+x_dist_2+x_dist_3)
      max_xx = maximum(x_dist_1+x_dist_2+x_dist_3);
      max_dd = maximum(d_c1s_Ek[:,2i]-reverse(z_baseline[(i-1)*NN+1:i*NN]));
      plot(d_c1s_Ek[:,(2i)-1],(max_xx/max_dd)*(d_c1s_Ek[:,2i]-reverse(z_baseline[(i-1)*NN+1:i*NN])))
   end
end


##
## cross section spread
##
Y1 = [reverse(d_c1s_Ek[:,2])-reverse(z_baseline_1); 1.0; zeros(Cdouble,Ny-2)];
Y2 = [reverse(d_c1s_Ek[:,4])-reverse(z_baseline_2); 1.0; zeros(Cdouble,Ny-2)];
Y3 = [reverse(d_c1s_Ek[:,6])-reverse(z_baseline_3); 1.0; zeros(Cdouble,Ny-2)];
Y4 = [reverse(d_c1s_Ek[:,8])-reverse(z_baseline_4); 1.0; zeros(Cdouble,Ny-2)];
Y5 = [reverse(d_c1s_Ek[:,10])-reverse(z_baseline_5); 1.0; zeros(Cdouble,Ny-2)];


# using XPSinv
# D_2nd = D2nd(N);
# Ny
dKe1 = d_c1s_Ek[1,1]-d_c1s_Ek[2,1];
dKe2 = d_c1s_Ek[1,3]-d_c1s_Ek[2,3];
dKe3 = d_c1s_Ek[1,5]-d_c1s_Ek[2,5];
dKe4 = d_c1s_Ek[1,7]-d_c1s_Ek[2,7];
dKe5 = d_c1s_Ek[1,9]-d_c1s_Ek[2,9];
R1 = dKe1*sum(d_c1s_Ek[:,2]-z_baseline_1);
R2 = dKe2*sum(d_c1s_Ek[:,4]-z_baseline_2);
R3 = dKe3*sum(d_c1s_Ek[:,6]-z_baseline_3);
R4 = dKe4*sum(d_c1s_Ek[:,8]-z_baseline_4);
R5 = dKe5*sum(d_c1s_Ek[:,10]-z_baseline_5);

Rreg1 = [R1*Matrix{Cdouble}(I,Ny,Ny); dKe1*ones(Cdouble,Ny)'; D_2nd];
Rreg2 = [R2*Matrix{Cdouble}(I,Ny,Ny); dKe2*ones(Cdouble,Ny)'; D_2nd];
Rreg3 = [R3*Matrix{Cdouble}(I,Ny,Ny); dKe3*ones(Cdouble,Ny)'; D_2nd];
Rreg4 = [R4*Matrix{Cdouble}(I,Ny,Ny); dKe4*ones(Cdouble,Ny)'; D_2nd];
Rreg5 = [R5*Matrix{Cdouble}(I,Ny,Ny); dKe5*ones(Cdouble,Ny)'; D_2nd];

ΓI1 = zeros(Cdouble,Ny+1+(Ny-2),Ny+1+(Ny-2));
ΓI2 = zeros(Cdouble,Ny+1+(Ny-2),Ny+1+(Ny-2));
ΓI3 = zeros(Cdouble,Ny+1+(Ny-2),Ny+1+(Ny-2));
ΓI4 = zeros(Cdouble,Ny+1+(Ny-2),Ny+1+(Ny-2));
ΓI5 = zeros(Cdouble,Ny+1+(Ny-2),Ny+1+(Ny-2));

σ_I1 = 100.0;
σ_I2 = 400.0;
σ_I3 = 500.0;
σ_I4 = 10.0;
σ_I5 = 20.0;

ΓI1[1:Ny,1:Ny] = σ_I1^2*Matrix{Cdouble}(I,Ny,Ny); # depends on the noise level
ΓI1[Ny+1,Ny+1] = (0.5/(d_c1s_Ek[1,1]-d_c1s_Ek[end,1]))^2; # 0.01^2
ΓI1[Ny+2:end,Ny+2:end] = 0.001^2*Matrix{Cdouble}(I,Ny-2,Ny-2);

ΓI2[1:Ny,1:Ny] = σ_I2^2*Matrix{Cdouble}(I,Ny,Ny); # depends on the noise level
# ΓI2[1:Ny,1:Ny] = 100.0^2*Matrix{Cdouble}(I,Ny,Ny);
ΓI2[Ny+1,Ny+1] = (0.5/(d_c1s_Ek[1,3]-d_c1s_Ek[end,3]))^2; # 0.01^2
ΓI2[Ny+2:end,Ny+2:end] = 0.001^2*Matrix{Cdouble}(I,Ny-2,Ny-2);

ΓI3[1:Ny,1:Ny] = σ_I3^2*Matrix{Cdouble}(I,Ny,Ny); # depends on the noise level
# ΓI3[1:Ny,1:Ny] = 100.0^2*Matrix{Cdouble}(I,Ny,Ny);
ΓI3[Ny+1,Ny+1] = (0.5/(d_c1s_Ek[1,5]-d_c1s_Ek[end,5]))^2; # 0.01^2
ΓI3[Ny+2:end,Ny+2:end] = 0.001^2*Matrix{Cdouble}(I,Ny-2,Ny-2);

ΓI4[1:Ny,1:Ny] = σ_I4^2*Matrix{Cdouble}(I,Ny,Ny); # depends on the noise level
# ΓI4[1:Ny,1:Ny] = 100.0^2*Matrix{Cdouble}(I,Ny,Ny);
ΓI4[Ny+1,Ny+1] = (0.5/(d_c1s_Ek[1,7]-d_c1s_Ek[end,7]))^2; # 0.01^2
ΓI4[Ny+2:end,Ny+2:end] = 0.001^2*Matrix{Cdouble}(I,Ny-2,Ny-2);

ΓI5[1:Ny,1:Ny] = σ_I5^2*Matrix{Cdouble}(I,Ny,Ny); # depends on the noise level
# ΓI5[1:Ny,1:Ny] = 100.0^2*Matrix{Cdouble}(I,Ny,Ny);
ΓI5[Ny+1,Ny+1] = (0.5/(d_c1s_Ek[1,9]-d_c1s_Ek[end,9]))^2; # 0.01^2
ΓI5[Ny+2:end,Ny+2:end] = 0.001^2*Matrix{Cdouble}(I,Ny-2,Ny-2);

ΓI1_inv = inv(ΓI1);
W1 = Rreg1'*ΓI1_inv*Rreg1;
W1_inv = inv(W1);

ΓI2_inv = inv(ΓI2);
W2 = Rreg2'*ΓI2_inv*Rreg2;
W2_inv = inv(W2);

ΓI3_inv = inv(ΓI3);
W3 = Rreg3'*ΓI3_inv*Rreg3;
W3_inv = inv(W3);

ΓI4_inv = inv(ΓI4);
W4 = Rreg4'*ΓI4_inv*Rreg4;
W4_inv = inv(W4);

ΓI5_inv = inv(ΓI5);
W5 = Rreg5'*ΓI5_inv*Rreg5;
W5_inv = inv(W5);


σ_CS1_reg = W1_inv*Rreg1'*Y1;
σ_CS2_reg = W2_inv*Rreg2'*Y2;
σ_CS3_reg = W3_inv*Rreg3'*Y3;
σ_CS4_reg = W4_inv*Rreg4'*Y4;
σ_CS5_reg = W5_inv*Rreg5'*Y5;

# σ_CS1_reg = reverse(d_c1s_Ek[:,2]-z_baseline_1)/R1;
figure(); plot(σ_CS1_reg); plot(σ_CS2_reg); plot(σ_CS3_reg); plot(σ_CS4_reg); plot(σ_CS5_reg);
# figure(); plot(reverse(d_c1s_Ek[:,2])-reverse(z_baseline_1))
# figure(); plot(Rreg1*σ_CS1_reg)

dKe1*sum(σ_CS1_reg)
dKe2*sum(σ_CS2_reg)
dKe3*sum(σ_CS3_reg)
dKe4*sum(σ_CS4_reg)
dKe5*sum(σ_CS5_reg)

using NewtonMethod
println("BFGS B")
function F1(x_::Array{Cdouble,1})
    0.5*(Y1-Rreg1*x_)'*ΓI1_inv*(Y1-Rreg1*x_)
end
function F2(x_::Array{Cdouble,1})
    0.5*(Y2-Rreg2*x_)'*ΓI2_inv*(Y2-Rreg2*x_)
end
function F3(x_::Array{Cdouble,1})
    0.5*(Y3-Rreg3*x_)'*ΓI3_inv*(Y3-Rreg3*x_)
end
function F4(x_::Array{Cdouble,1})
    0.5*(Y4-Rreg4*x_)'*ΓI4_inv*(Y4-Rreg4*x_)
end
function F5(x_::Array{Cdouble,1})
    0.5*(Y5-Rreg5*x_)'*ΓI5_inv*(Y5-Rreg5*x_)
end

function Fgrad1(x_::Array{Cdouble,1})
    -Rreg1'*ΓI1_inv*(Y1-Rreg1*x_)
end
function Fgrad2(x_::Array{Cdouble,1})
    -Rreg2'*ΓI2_inv*(Y2-Rreg2*x_)
end
function Fgrad3(x_::Array{Cdouble,1})
    -Rreg3'*ΓI3_inv*(Y3-Rreg3*x_)
end
function Fgrad4(x_::Array{Cdouble,1})
    -Rreg4'*ΓI4_inv*(Y4-Rreg4*x_)
end
function Fgrad5(x_::Array{Cdouble,1})
    -Rreg5'*ΓI5_inv*(Y5-Rreg5*x_)
end

"""
I_nbl:   photo-electric signal (corrected for the baseline)
Kes:     kinetic energy discretization points (regularly sub division)
σ_I:     noise level in the measurement (standard deviation)
Nbfgs:   number of iterations in the bounded BFGS loop
Nsearch: maximum number of iteration for the line search
"""
function cross_section_spread_function(I_nbl::Array{Cdouble,1},Kes::Array{Cdouble,1},σ_I::Cdouble;Nbfgs::Int64=1000,Nsearch::Int64=10)
   Ny = length(I_nbl)
   if (length(Kes)!=Ny) throw("not the right amount of kinetic energies") end
   dKe = Kes[2] - Kes[1];
   R = dKe*sum(I_nbl);
   D_2nd = diagm(Ny-2,Ny,1 => 2ones(Cdouble,Ny-2), 0 => -ones(Cdouble,Ny-2) ,2 => -ones(Cdouble,Ny-2));

   # augmented measurement
   Y = [I_nbl; 1.0; zeros(Cdouble,Ny-2)];

   # augmented measurement operator (augmented => a priori in the operator)
   Rreg = [R*Matrix{Cdouble}(I,Ny,Ny); dKe*ones(Cdouble,Ny)'; D_2nd];

   # covariance matrix of the augmented measurements
   ΓI = zeros(Cdouble,Ny+1+(Ny-2),Ny+1+(Ny-2));
   ΓI[1:Ny,1:Ny] = σ_I^2*Matrix{Cdouble}(I,Ny,Ny);                # noise level in the measurment
   ΓI[Ny+1,Ny+1] = (0.5/(Kes[1]-Kes[end]))^2;                     # accuracy of the numerical integration
   ΓI[Ny+2:end,Ny+2:end] = 0.001^2*Matrix{Cdouble}(I,Ny-2,Ny-2);  # variance of the second order difference of the CS spread function
   ΓI_inv = inv(ΓI);
   W_inv  = inv(Rreg'*ΓI_inv*Rreg);

   # optimization problem
   function F(x::Array{Cdouble,1})
       0.5*(Y-Rreg*x)'*ΓI_inv*(Y-Rreg*x)
   end

   function Fgrad(x::Array{Cdouble,1})
       -Rreg'*ΓI_inv*(Y-Rreg*x)
   end

   X0 = zeros(Cdouble,Ny);
   p0 = zeros(Cdouble,Ny);            # first descent direction
   alpha_min = -4.0                   # smallest value of the length of the step
   alpha_max = 4.0                    # smallest value of the length of the step 2000000.0
   mu = 0.4                           # <0.5 parameter for the line search algorithm
   # H0 = Matrix{Cdouble}(I,Ny,Ny)/R1;  # initial inverse Hessian matrix
   H0 = W_inv
   # define the constraints
   lx = zeros(Cdouble,Ny);            # lower bounds
   ux = Inf*ones(Cdouble,Ny);         # upper bounds
   # Xend,Hend,Xpath,Nlast = BFGSB(X0,H0,Nbfgs,alpha_min,alpha_max,mu,lx,ux,F,Fgrad,Nsearch);
   BFGSB(X0,H0,Nbfgs,alpha_min,alpha_max,mu,lx,ux,F,Fgrad,Nsearch)
end
X0 = zeros(Cdouble,Ny); # 0.5σ_CS1_reg #
p0 = zeros(Cdouble,Ny); # Fgrad(0.5σ_CS1_reg); # ones(Cdouble,Ny);             # first descent direction
alpha_min = -4.0          # smallest value of the length of the step
alpha_max = 4.0          # smallest value of the length of the step 2000000.0
mu = 0.4                 # <0.5 parameter for the line search algorithm
Nbfgs = 1000             # number of iterations
Nsearch = 10             # maximum number of iteration for the line search
H01 = Matrix{Cdouble}(I,Ny,Ny)/R1; # W1_inv;                   # initial Hessian matrix #0.000001*
H02 = Matrix{Cdouble}(I,Ny,Ny)/R2;
H03 = Matrix{Cdouble}(I,Ny,Ny)/R3;
H04 = Matrix{Cdouble}(I,Ny,Ny)/R4;
H05 = Matrix{Cdouble}(I,Ny,Ny)/R5;

# define the constraints
lx = zeros(Cdouble,Ny);
ux = Inf*ones(Cdouble,Ny);
# Xend1,Hend1,Xpath1,Nlast1 = BFGSB(X0,H01,Nbfgs,alpha_min,alpha_max,mu,lx,ux,F,Fgrad1,Nsearch);
Xend1,Hend1,Xpath1,Nlast1 = cross_section_spread_function(reverse(d_c1s_Ek[:,2])-reverse(z_baseline_1),reverse(d_c1s_Ek[:,1]),σ_I1;Nbfgs=1000,Nsearch=10)
# Xend2,Hend2,Xpath2,Nlast2 = BFGSB(X0,H02,Nbfgs,alpha_min,alpha_max,mu,lx,ux,F2,Fgrad2,Nsearch);
Xend2,Hend2,Xpath2,Nlast2 = cross_section_spread_function(reverse(d_c1s_Ek[:,4])-reverse(z_baseline_2),reverse(d_c1s_Ek[:,3]),σ_I2;Nbfgs=1000,Nsearch=10)
# Xend3,Hend3,Xpath3,Nlast3 = BFGSB(X0,H03,Nbfgs,alpha_min,alpha_max,mu,lx,ux,F3,Fgrad3,Nsearch);
Xend3,Hend3,Xpath3,Nlast3 = cross_section_spread_function(reverse(d_c1s_Ek[:,6])-reverse(z_baseline_3),reverse(d_c1s_Ek[:,5]),σ_I3;Nbfgs=1000,Nsearch=10)
# Xend4,Hend4,Xpath4,Nlast4 = BFGSB(X0,H04,Nbfgs,alpha_min,alpha_max,mu,lx,ux,F4,Fgrad4,Nsearch);
Xend4,Hend4,Xpath4,Nlast4 = cross_section_spread_function(reverse(d_c1s_Ek[:,8])-reverse(z_baseline_4),reverse(d_c1s_Ek[:,7]),σ_I4;Nbfgs=1000,Nsearch=10)
# Xend5,Hend5,Xpath5,Nlast5 = BFGSB(X0,H05,Nbfgs,alpha_min,alpha_max,mu,lx,ux,F5,Fgrad5,Nsearch);
Xend5,Hend5,Xpath5,Nlast5 = cross_section_spread_function(reverse(d_c1s_Ek[:,10])-reverse(z_baseline_5),reverse(d_c1s_Ek[:,9]),σ_I5;Nbfgs=1000,Nsearch=10)
dKe1*sum(Xend1)
dKe2*sum(Xend2)
dKe3*sum(Xend3)
dKe4*sum(Xend4)
dKe5*sum(Xend5)
figure(); plot(Xend1); plot((reverse(d_c1s_Ek[:,2])-reverse(z_baseline_1))/R1)
figure(); plot(Xend2); plot((reverse(d_c1s_Ek[:,4])-reverse(z_baseline_2))/R2)
figure(); plot(Xend3); plot((reverse(d_c1s_Ek[:,6])-reverse(z_baseline_3))/R3)
figure(); plot(Xend4); plot((reverse(d_c1s_Ek[:,8])-reverse(z_baseline_4))/R4)
figure(); plot(Xend5); plot((reverse(d_c1s_Ek[:,10])-reverse(z_baseline_5))/R5)

figure();
plot(reverse(d_c1s_Ek[:,1]),Xend1)
plot(reverse(d_c1s_Ek[:,3]),Xend2)
plot(reverse(d_c1s_Ek[:,5]),Xend3)
plot(reverse(d_c1s_Ek[:,7]),Xend4)
plot(reverse(d_c1s_Ek[:,9]),Xend5)


##
## load the computed integrals
##
d_Na = CSV.File(string("data/","na_prop_prop_2_output(1).csv"); delim=",", header=true) |> DataFrame


## number of integration points for the Simpsons rule
Nz0 = 50;
## ρ_tot_int: total concentration. It should be mainly just water concentration
σ_z0 = 0.5; # [5 Å] width of the transition region
z00 = 0.5;  # half height depth

ħν_exp = 1.0*d_Na[1:3:end,5];
Fν_exp = 1.0*d_Na[1:3:end,4]; # 1.0e-9
#WARNING: rectifying mirror current, maybe the gain to transform the current to photon flux is not the same for all photon energy!!!!!
Fν_exp[4] = Fν_exp[4]/20.0;
Fν_exp[5] = Fν_exp[5]/45.0;
T_exp  = d_Na[1:3:end,9];
μKe_exp = d_Na[1:3:end,8];
Be_exp = [d_c1s_Ek[:,1]'; d_c1s_Ek[:,3]'; d_c1s_Ek[:,5]'; d_c1s_Ek[:,7]'; d_c1s_Ek[:,9]'];
σν_exp_1 = zeros(Cdouble,5,length(d_c1s_Ek[:,1]));
σν_exp_2 = zeros(Cdouble,5,length(d_c1s_Ek[:,1]));
σν_exp_3 = zeros(Cdouble,5,length(d_c1s_Ek[:,1]));


##
##TODO: use one peak by one peak to create models and then data. Later, create bad models in the last two channels
##
for i in 1:5
   x_dist_1 = (τt[i,1]/σt[i,1])*exp.(-0.5*((d_c1s_Ek[:,(2i)-1].-μt[i,1])/σt[i,1]).^2);
   x_dist_2 = (τt[i,2]/σt[i,2])*exp.(-0.5*((d_c1s_Ek[:,(2i)-1].-μt[i,2])/σt[i,2]).^2);
   x_dist_3 = (τt[i,3]/σt[i,3])*exp.(-0.5*((d_c1s_Ek[:,(2i)-1].-μt[i,3])/σt[i,3]).^2);
   max_xx = maximum(x_dist_1+x_dist_2+x_dist_3);
   σν_exp_1[i,:] = (σ_cs_orb(ħν_exp[i],"C1s")/max_xx)*x_dist_1
   σν_exp_2[i,:] = (σ_cs_orb(ħν_exp[i],"C1s")/max_xx)*x_dist_2
   σν_exp_3[i,:] = (σ_cs_orb(ħν_exp[i],"C1s")/max_xx)*x_dist_3
end

Fν_exp   .= 1.0;
T_exp    .= 1.0;
σν_exp_1 = σν_exp_1./σ_cs_orb.(ħν_exp,"C1s");
σν_exp_2 = σν_exp_2./σ_cs_orb.(ħν_exp,"C1s");
σν_exp_3 = σν_exp_3./σ_cs_orb.(ħν_exp,"C1s");

wsXPS_1 = XPSsetup(ħν_exp,Fν_exp,μKe_exp,T_exp,Be_exp,σν_exp_1;α_exp=1.0);
wsXPS_2 = XPSsetup(ħν_exp,Fν_exp,μKe_exp,T_exp,Be_exp,σν_exp_2;α_exp=1.0);
wsXPS_3 = XPSsetup(ħν_exp,Fν_exp,μKe_exp,T_exp,Be_exp,σν_exp_3;α_exp=1.0);


# depth discretization
N = 50;
Z_max = 10.0;
Zi = collect(range(0.0,Z_max,length=N));
H_1 = Ψ_lin_peaks(Zi,wsXPS_1;Nz=Nz0,σ_z=σ_z0,z0=z00,κ_cs=0.0,κ_eal=0.0);
H_2 = Ψ_lin_peaks(Zi,wsXPS_2;Nz=Nz0,σ_z=σ_z0,z0=z00,κ_cs=0.0,κ_eal=0.0);
H_3 = Ψ_lin_peaks(Zi,wsXPS_3;Nz=Nz0,σ_z=σ_z0,z0=z00,κ_cs=0.0,κ_eal=0.0);
# figure(); imshow(H_1); colorbar()
# figure(); imshow(H_2); colorbar()
# figure(); imshow(H_3); colorbar()
H_mean_1,H_std_1 = Ψ_lin_peaks_mean_and_std(Zi,wsXPS_1;Nz=Nz0,κ_cs=0.05,κ_eal=0.05,σ_z=σ_z0,z0=z00);
H_mean_2,H_std_2 = Ψ_lin_peaks_mean_and_std(Zi,wsXPS_2;Nz=Nz0,κ_cs=0.05,κ_eal=0.05,σ_z=σ_z0,z0=z00);
H_mean_3,H_std_3 = Ψ_lin_peaks_mean_and_std(Zi,wsXPS_3;Nz=Nz0,κ_cs=0.05,κ_eal=0.05,σ_z=σ_z0,z0=z00);

##
## profiles
##

ρA_1 = logistic.(Zi.-2.0,0.0,1.0,2.0);
ρA_2 = logistic.(Zi.-2.0,0.0,1.0,2.0) .+ 2.0exp.(-(Zi.-1.0).^2. /(2.0*0.25^2));
ρA_3 = logistic.(Zi.-2.0,0.0,1.0,2.0) .+ exp.(-(Zi.-1.5).^2. /(2.0*0.5^2));
ρA_4 = exp.(-(Zi.-2.5).^2. /(2.0*0.5^2));


##
## generate spectra
##
IA_1_clean_1 = H_1*ρA_1;
IA_2_clean_1 = H_1*ρA_2;
IA_3_clean_1 = H_1*ρA_3;
IA_4_clean_1 = H_1*ρA_4;

σ_noise = 0.0
IA_1_1 = IA_1_clean_1 + σ_noise*randn(wsXPS_1.Nke*wsXPS_1.Nbe);
IA_2_1 = IA_2_clean_1 + σ_noise*randn(wsXPS_1.Nke*wsXPS_1.Nbe);
IA_3_1 = IA_3_clean_1 + σ_noise*randn(wsXPS_1.Nke*wsXPS_1.Nbe);
IA_4_1 = IA_4_clean_1 + σ_noise*randn(wsXPS_1.Nke*wsXPS_1.Nbe);

figure();
for i in 1:wsXPS_1.Nν
   scatter(wsXPS_1.Be[i,:],IA_1_1[(i-1)*wsXPS_1.Nbe+1:i*wsXPS_1.Nbe]) # ,color="tab:blue")
   plot(wsXPS_1.Be[i,:],IA_1_clean_1[(i-1)*wsXPS_1.Nbe+1:i*wsXPS_1.Nbe]) # ,color="tab:blue")
end
xlabel("B\$_e\$ [a.u.]")
ylabel("PE signal [a.u.]")
xlim(minimum(wsXPS_1.Be),maximum(wsXPS_1.Be))
title("simulated spectra")
ylim(0.0,1.1maximum(IA_1_1))


figure();
for i in 1:wsXPS_1.Nν
   scatter(wsXPS_1.Be[i,:],IA_2_1[(i-1)*wsXPS_1.Nbe+1:i*wsXPS_1.Nbe]) # ,color="tab:blue")
   plot(wsXPS_1.Be[i,:],IA_2_clean_1[(i-1)*wsXPS_1.Nbe+1:i*wsXPS_1.Nbe]) # ,color="tab:blue")
end
xlabel("B\$_e\$ [a.u.]")
ylabel("PE signal [a.u.]")
xlim(minimum(wsXPS_1.Be),maximum(wsXPS_1.Be))
title("simulated spectra")
ylim(0.0,1.1maximum(IA_2_1))


figure();
for i in 1:wsXPS_1.Nν
   scatter(wsXPS_1.Be[i,:],IA_3_1[(i-1)*wsXPS_1.Nbe+1:i*wsXPS_1.Nbe]) # ,color="tab:blue")
   plot(wsXPS_1.Be[i,:],IA_3_clean_1[(i-1)*wsXPS_1.Nbe+1:i*wsXPS_1.Nbe]) # ,color="tab:blue")
end
xlabel("B\$_e\$ [a.u.]")
ylabel("PE signal [a.u.]")
xlim(minimum(wsXPS_1.Be),maximum(wsXPS_1.Be))
title("simulated spectra")
ylim(0.0,1.1maximum(IA_3_1))


figure();
for i in 1:wsXPS_1.Nν
   scatter(wsXPS_1.Be[i,:],IA_4_1[(i-1)*wsXPS_1.Nbe+1:i*wsXPS_1.Nbe]) # ,color="tab:blue")
   plot(wsXPS_1.Be[i,:],IA_4_clean_1[(i-1)*wsXPS_1.Nbe+1:i*wsXPS_1.Nbe]) # ,color="tab:blue")
end
xlabel("B\$_e\$ [a.u.]")
ylabel("PE signal [a.u.]")
xlim(minimum(wsXPS_1.Be),maximum(wsXPS_1.Be))
title("simulated spectra")
ylim(0.0,1.1maximum(IA_4_1))

if SAVE_MODEL
   # sample the relative error in cross section and eal to generate measurement model with modelling error
   model_folder = "./data/lin/real/low_res_50/more/";
   κ_css  = collect(range(-0.05,0.05,length=11));
   κ_eals = collect(range(-0.05,0.05,length=11));
   i = 0
   t_elapsed = @elapsed for (κ1,κ2) in Iterators.product(κ_css,κ_eals)
      global i = i + 1;
      # compute the measurment model with some error
      # H,λeal,σ_cs = Ψ_lin_peaks(Zi,Kes,μ_σ_ke;Nz=Nz0,σ_z=σ_z0,z0=z00,κ_cs=κ1,κ_eal=κ2);
      Hκ = Ψ_lin_peaks(Zi,wsXPS;Nz=Nz0,σ_z=σ_z0,z0=z00,κ_cs=κ1,κ_eal=κ2);
      # save the model, the intended and true crosss section and eal
      s = @sprintf "%s/%i/" model_folder i
      mkpath(s)
      CSV.write(string(s,"kinetic_energy.csv"),DataFrame(wsXPS.Ke',:auto); header=false);
      CSV.write(string(s,"depth.csv"),DataFrame(Zi',:auto); header=false);
      CSV.write(string(s,"H.csv"),DataFrame(Hκ,:auto); header=false);
   end
end

if false
   # rescale H
   # y_all = [d_c1s_Ek[:,2]; d_c1s_Ek[:,4]; d_c1s_Ek[:,6]; d_c1s_Ek[:,8];d_c1s_Ek[:,10]]
   y_all = [d_c1s_Ek[:,2].-z_baseline_1; d_c1s_Ek[:,4].-z_baseline_2; d_c1s_Ek[:,6].-z_baseline_3; d_c1s_Ek[:,8].-z_baseline_4;d_c1s_Ek[:,10].-z_baseline_5]
   H = (maximum(y_all[4wsXPS.Nbe+1:5wsXPS.Nbe])/maximum(H[4wsXPS.Nbe+1:5wsXPS.Nbe,2]))*H

   F = svd(H, full=true);
   F1 = svd(H[1:wsXPS.Nbe,:], full=true);
   F2 = svd(H[wsXPS.Nbe+1:2wsXPS.Nbe,:], full=true);
   F3 = svd(H[2wsXPS.Nbe+1:3wsXPS.Nbe,:], full=true);
   F4 = svd(H[3wsXPS.Nbe+1:4wsXPS.Nbe,:], full=true);
   F5 = svd(H[4wsXPS.Nbe+1:5wsXPS.Nbe,:], full=true);

   abs_ub_1            = abs.(F.U'*IA_1);
   abs_ub_2            = abs.(F.U'*IA_2);
   abs_ub_3            = abs.(F.U'*IA_3);
   min_ub,max_ub = extrema([abs_ub_1; abs_ub_2; abs_ub_3])
   l_scat_4_ub_1        = scatter(collect(1:wsXPS.Nbe*wsXPS.Nν),abs.(abs_ub_1),color="tab:blue")
   l_scat_4_ub_2        = scatter(collect(1:wsXPS.Nbe*wsXPS.Nν),abs.(abs_ub_2),color="tab:orange")
   l_scat_4_ub_3        = scatter(collect(1:wsXPS.Nbe*wsXPS.Nν),abs.(abs_ub_3),color="tab:pink")
   ylim(0.5min_ub,5.0max_ub)
   yscale("log")
   xlabel("SVD left index",fontsize=12)
   ylabel("SVD left coefficient",fontsize=12)
   legend([l_scat_4_ub_1,l_scat_4_ub_2,l_scat_4_ub_3],["data 1","data 2","darta 3"])


   figure()
   abs_ub            = abs.(F.U'*y_all);
   min_ub,max_ub = extrema(abs_ub)
   l_scat_4_ub        = scatter(collect(1:wsXPS.Nbe*wsXPS.Nν),abs.(abs_ub),color="tab:blue")
   ylim(0.5min_ub,5.0max_ub)
   yscale("log")
   xlabel("SVD left index",fontsize=12)
   ylabel("SVD left coefficient",fontsize=12)
   legend([l_scat_4_ub],["data"])

   figure()
   abs_ub_1 = abs.(F1.U'*(d_c1s_Ek[:,2]-z_baseline_1))
   abs_ub_2 = abs.(F2.U'*(d_c1s_Ek[:,4]-z_baseline_2))
   abs_ub_3 = abs.(F3.U'*(d_c1s_Ek[:,6]-z_baseline_3))
   abs_ub_4 = abs.(F4.U'*(d_c1s_Ek[:,8]-z_baseline_4))
   abs_ub_5 = abs.(F5.U'*(d_c1s_Ek[:,10]-z_baseline_5))
   σ_1 = mean(abs_ub_1[50:end])
   σ_2 = mean(abs_ub_2[50:end])
   σ_3 = mean(abs_ub_3[50:end])
   σ_4 = mean(abs_ub_4[50:end])
   σ_5 = mean(abs_ub_5[50:end])
   min_ub,max_ub = extrema([abs_ub_1; abs_ub_2; abs_ub_3; abs_ub_4; abs_ub_5;])
   l_scat_ub_1        = scatter(collect(1:wsXPS.Nbe),abs.(abs_ub_1),color="tab:blue")
   plot([1;wsXPS.Nbe], [σ_1;σ_1],color="tab:blue")
   l_scat_ub_2        = scatter(collect(wsXPS.Nbe+1:2wsXPS.Nbe),abs.(abs_ub_2),color="tab:orange")
   plot([wsXPS.Nbe+1;2wsXPS.Nbe], [σ_2;σ_2],color="tab:orange")
   l_scat_ub_3        = scatter(collect(2wsXPS.Nbe+1:3wsXPS.Nbe),abs.(abs_ub_3),color="tab:green")
   plot([2wsXPS.Nbe+1;3wsXPS.Nbe], [σ_3;σ_3],color="tab:green")
   l_scat_ub_4        = scatter(collect(3wsXPS.Nbe+1:4wsXPS.Nbe),abs.(abs_ub_4),color="tab:red")
   plot([3wsXPS.Nbe+1;4wsXPS.Nbe], [σ_4;σ_4],color="tab:red")
   l_scat_ub_5        = scatter(collect(4wsXPS.Nbe+1:5wsXPS.Nbe),abs.(abs_ub_5),color="tab:cyan")
   plot([4wsXPS.Nbe+1;5wsXPS.Nbe], [σ_5;σ_5],color="tab:cyan")
   ylim(0.5min_ub,5.0max_ub)
   yscale("log")
   xlabel("SVD left index",fontsize=12)
   ylabel("SVD left coefficient",fontsize=12)
   legend([l_scat_ub_1,l_scat_ub_2,l_scat_ub_3,l_scat_ub_4,l_scat_ub_5],["data 1"; "data 2"; "data 3"; "data 4"; "data 5"])






   figure()
   l_scat_1        = scatter(collect(1:wsXPS.Nbe),d_c1s_Ek[:,2]-z_baseline_1,color="tab:blue")
   l_scat_2        = scatter(collect(wsXPS.Nbe+1:2wsXPS.Nbe),d_c1s_Ek[:,4]-z_baseline_2,color="tab:orange")
   l_scat_3        = scatter(collect(2wsXPS.Nbe+1:3wsXPS.Nbe),d_c1s_Ek[:,6]-z_baseline_3,color="tab:green")
   l_scat_4        = scatter(collect(3wsXPS.Nbe+1:4wsXPS.Nbe),d_c1s_Ek[:,8]-z_baseline_4,color="tab:red")
   l_scat_5        = scatter(collect(4wsXPS.Nbe+1:5wsXPS.Nbe),d_c1s_Ek[:,10]-z_baseline_5,color="tab:cyan")
   xlabel("SVD left index",fontsize=12)
   ylabel("PE signal [a.u.]",fontsize=12)
   legend([l_scat_1,l_scat_2,l_scat_3,l_scat_4,l_scat_5],["data 1"; "data 2"; "data 3"; "data 4"; "data 5"])



   figure()
   plot(wsXPS.Ke[4wsXPS.Nbe+1:5wsXPS.Nbe],y_all[4wsXPS.Nbe+1:5wsXPS.Nbe])
   plot(wsXPS.Ke[4wsXPS.Nbe+1:5wsXPS.Nbe],(maximum(y_all[4wsXPS.Nbe+1:5wsXPS.Nbe])/maximum(H[4wsXPS.Nbe+1:5wsXPS.Nbe,2]))*H[4wsXPS.Nbe+1:5wsXPS.Nbe,2])


   figure()
   plot(wsXPS.Ke[wsXPS.Nbe+1:2wsXPS.Nbe],y_all[wsXPS.Nbe+1:2wsXPS.Nbe])
   plot(wsXPS.Ke[wsXPS.Nbe+1:2wsXPS.Nbe],(maximum(y_all[wsXPS.Nbe+1:2wsXPS.Nbe])/maximum(H[wsXPS.Nbe+1:2wsXPS.Nbe,2]))*H[wsXPS.Nbe+1:2wsXPS.Nbe,2])


   minH,maxH = extrema(H)
   fig1, ax1, pcm1 = imshowData(1254,Zi,wsXPS.Ke[wsXPS.Nbe+1:2wsXPS.Nbe],H[wsXPS.Nbe+1:2wsXPS.Nbe,:],_norm=:Normalize,_vmin=minH,_vmax=maxH,_edgecolors="face",_sub=221)
   xlabel("depth \$z\$",fontsize=14) # [nm]
   ylabel("Kinetic energy \$K_e\$",fontsize=14)
   xlim(Zi[1],Zi[end])
   ylim(wsXPS.Ke[wsXPS.Nbe+1],wsXPS.Ke[2wsXPS.Nbe])
   rc("ytick",color="white")
   cax1 = fig1.add_axes([0.08+0.33, .75, 0.02, 0.2])
   cb1 = fig1.colorbar(pcm1, orientation="vertical", cax=cax1, shrink=0.6)
   cb1.set_label("measurement operator [a.u.]", color="white", fontsize=10)
   cb1.ax.yaxis.set_tick_params(color="white")
   cb1.outline.set_edgecolor("white")
   rc("ytick",color="black")


   minU,maxU = extrema(F.U)
   minU,maxU = extrema(F.U[wsXPS.Nbe+1:2wsXPS.Nbe,:])
   fig1, ax2, pcm2 = imshowData(1254,collect(1.0:wsXPS.Nbe*wsXPS.Nν),wsXPS.Ke[wsXPS.Nbe+1:2wsXPS.Nbe],F.U[wsXPS.Nbe+1:2wsXPS.Nbe,:],_norm=:Normalize,_vmin=minU,_vmax=maxU,_edgecolors="face",_sub=222)
   xlabel("basis index [data space]",fontsize=14) # [nm]
   ylabel("Kinetic energy \$K_e\$",fontsize=14)
   xlim(1,length(wsXPS.Ke))
   ylim(wsXPS.Ke[wsXPS.Nbe+1],wsXPS.Ke[2wsXPS.Nbe])
   rc("ytick",color="white")
   cax2 = fig1.add_axes([0.57+0.31, .77, 0.02, 0.2])
   cb2 = fig1.colorbar(pcm2, orientation="vertical", cax=cax2, shrink=0.6)
   cb2.set_label("left vector entries [a.u.]", color="white", fontsize=10)
   cb2.ax.yaxis.set_tick_params(color="white")
   cb2.outline.set_edgecolor("white")
   rc("ytick",color="black")


   minV,maxV = extrema(F.Vt)
   fig1, ax3, pcm3 = imshowData(1254,Zi,collect(1.0:N),F.Vt,_norm=:Normalize,_vmin=minV,_vmax=maxV,_edgecolors="face",_sub=223)
   xlabel("depth \$z\$",fontsize=14) # [nm]
   ylabel("basis index [depth space]",fontsize=14)
   xlim(Zi[1],Zi[end])
   ylim(1,N)
   rc("ytick",color="white")
   cax3 = fig1.add_axes([0.08, .27, 0.02, 0.2])
   cb3 = fig1.colorbar(pcm3, orientation="vertical", cax=cax3, shrink=0.6)
   cb3.set_label("right vector entries [a.u.]", color="white", fontsize=10)
   cb3.ax.yaxis.set_tick_params(color="white")
   cb3.outline.set_edgecolor("white")
   rc("ytick",color="black")


   ax4 = subplot(224)
   semilogy(F.S)
   xlabel("basis index [depth space]",fontsize=14)
   ylabel("singular values [a.u.]",fontsize=14)


   fig1.set_figwidth(10.0) # (5.6)
   fig1.set_figheight(9.0)
   tight_layout(pad=1.0, w_pad=0.2, h_pad=0.2)

   ax3.annotate("a)", xy=(3, 1),  xycoords="data", xytext=(-0.12, 2.12), textcoords="axes fraction", color="black",fontsize=14)
   ax3.annotate("b)", xy=(3, 1),  xycoords="data", xytext=(1.07, 2.12), textcoords="axes fraction", color="black",fontsize=14)
   ax3.annotate("c)", xy=(3, 1),  xycoords="data", xytext=(-0.12, 1.0), textcoords="axes fraction", color="black",fontsize=14)
   ax3.annotate("d)", xy=(3, 1),  xycoords="data", xytext=(1.07, 1.0), textcoords="axes fraction", color="black",fontsize=14)
end







##
## try some inversion
##
if SOME_INV
   using XPSinv


   D_2nd = D2nd(N);
   d_2nd = 2;
   B = zeros(Cdouble,2,N); B[1,1] = 1.0; B[2,end] = 1.0
   γ = σ_1^2*ones(Cdouble,wsXPS.Nbe*wsXPS.Nke);
   γ[wsXPS.Nbe+1:2wsXPS.Nbe] .= σ_2^2;
   γ[2wsXPS.Nbe+1:3wsXPS.Nbe] .= σ_3^2;
   γ[3wsXPS.Nbe+1:4wsXPS.Nbe] .= σ_4^2;
   γ[4wsXPS.Nbe+1:5wsXPS.Nbe] .= σ_5^2;
   Γ_inv = diagm(1.0./γ);
   # z00  = 2.0;
   σ_z0 = 1.0;

   γ_D = logistic.(Zi.-3.0,0.2,0.0,2.0σ_z0)[2:end-1];
   # γ_D = 0.1logistic.(Zi.-3.0,0.2,0.0,2.0σ_z0)[2:end-1].+10.0;
   γ_D = (γ_D .+ 0.01maximum(γ_D)).^2;
   Γ_D_inv = diagm(1.0./γ_D);

   ρ0 = ρA_1[1] # 0.0;
   γ0 = 0.01^2 # 1.0^2 # 0.01^2 # 0.1^2
   ρB = ρA_1[end] # 1.0;
   γB = 0.01^2 # 0.01^2 # 0.1^2
   ρ_bar = [ρ0;ρB];
   γ_ρ = [γ0;γB];
   Γ_ρ_inv = diagm(1.0./γ_ρ);
   rankK = rank(H);



   ##
   ## draw samples from data distribution
   ##
   # μ_ρ,σ_ρ,ρ_samples,F,W_inv = data_sample_nso_un(1000,rankK,H,0.0H,σ_noise*ones(Cdouble,Nke),Γ_inv,D_2nd,Γ_D_inv,B,Γ_ρ_inv,IA_1_clean,ρ_bar)

   # ρ_clean,_,_ = iterative_nso(rankK,H,Γ_inv,D_2nd,Γ_D_inv,B,Γ_ρ_inv,y_all,ρ_bar);
   # ρ_clean,_,_ = iterative_nso(rankK,H,Γ_inv,D_2nd,Γ_D_inv,B,Γ_ρ_inv,H*ρA_1,ρ_bar);


   # ρ_clean,σ_ρ,ρ_samples,F,W_inv = data_sample_nso_un(100,rankK,H,0.01H,sqrt.(γ),Γ_inv,D_2nd,Γ_D_inv,B,Γ_ρ_inv,H*ρA_1,ρ_bar)
   Nmin = 50; # 20
   ρ_clean,σ_ρ,ρ_samples,F,W_inv = shuffle_data_sample_nso_un(100,wsXPS.Nke,rankK,H,1.0e-2H,γ,Γ_inv,D_2nd,Γ_D_inv,B,Γ_ρ_inv,H*ρA_1,ρ_bar;Nmin=Nmin);


   figure(); plot(Zi,ρ_clean)
   figure(); plot(Zi,ρ_samples')

   # ρ_clean .= 1.0
   # ρ_clean[:] = ρA_3/1.5
   #
   # ρ_clean[:] = 1.2ρA_1
   #
   #
   # ρ_clean[:] = ρA_3-0.9ρA_1

   figure()
   l_scat_1        = scatter(collect(1:wsXPS.Nbe),d_c1s_Ek[:,2]-z_baseline_1,color="tab:blue")
   l_plot_1,       = plot(collect(1:wsXPS.Nbe),H[1:wsXPS.Nbe,:]*ρ_clean,color="tab:blue")
   l_fill_1        = fill_between(collect(1:wsXPS.Nbe),H[1:wsXPS.Nbe,:]*ρ_clean-σ_1*ones(Cdouble,wsXPS.Nbe),H[1:wsXPS.Nbe,:]*ρ_clean+σ_1*ones(Cdouble,wsXPS.Nbe),alpha=0.5,color="tab:blue")
   l_scat_2        = scatter(collect(wsXPS.Nbe+1:2wsXPS.Nbe),d_c1s_Ek[:,4]-z_baseline_2,color="tab:orange")
   l_plot_2,       = plot(collect(wsXPS.Nbe+1:2wsXPS.Nbe),H[wsXPS.Nbe+1:2wsXPS.Nbe,:]*ρ_clean,color="tab:orange")
   l_fill_2        = fill_between(collect(wsXPS.Nbe+1:2wsXPS.Nbe),H[wsXPS.Nbe+1:2wsXPS.Nbe,:]*ρ_clean-σ_2*ones(Cdouble,wsXPS.Nbe),H[wsXPS.Nbe+1:2wsXPS.Nbe,:]*ρ_clean+σ_2*ones(Cdouble,wsXPS.Nbe),alpha=0.5,color="tab:orange")
   l_scat_3        = scatter(collect(2wsXPS.Nbe+1:3wsXPS.Nbe),d_c1s_Ek[:,6]-z_baseline_3,color="tab:green")
   l_plot_3,       = plot(collect(2wsXPS.Nbe+1:3wsXPS.Nbe),H[2wsXPS.Nbe+1:3wsXPS.Nbe,:]*ρ_clean,color="tab:green")
   l_fill_3        = fill_between(collect(2wsXPS.Nbe+1:3wsXPS.Nbe),H[2wsXPS.Nbe+1:3wsXPS.Nbe,:]*ρ_clean-σ_3*ones(Cdouble,wsXPS.Nbe),H[2wsXPS.Nbe+1:3wsXPS.Nbe,:]*ρ_clean+σ_3*ones(Cdouble,wsXPS.Nbe),alpha=0.5,color="tab:green")
   l_scat_4        = scatter(collect(3wsXPS.Nbe+1:4wsXPS.Nbe),d_c1s_Ek[:,8]-z_baseline_4,color="tab:red")
   l_plot_4,       = plot(collect(3wsXPS.Nbe+1:4wsXPS.Nbe),H[3wsXPS.Nbe+1:4wsXPS.Nbe,:]*ρ_clean,color="tab:red")
   l_fill_4        = fill_between(collect(3wsXPS.Nbe+1:4wsXPS.Nbe),H[3wsXPS.Nbe+1:4wsXPS.Nbe,:]*ρ_clean-σ_4*ones(Cdouble,wsXPS.Nbe),H[3wsXPS.Nbe+1:4wsXPS.Nbe,:]*ρ_clean+σ_4*ones(Cdouble,wsXPS.Nbe),alpha=0.5,color="tab:red")
   l_scat_5        = scatter(collect(4wsXPS.Nbe+1:5wsXPS.Nbe),d_c1s_Ek[:,10]-z_baseline_5,color="tab:cyan")
   l_plot_5,       = plot(collect(4wsXPS.Nbe+1:5wsXPS.Nbe),H[4wsXPS.Nbe+1:5wsXPS.Nbe,:]*ρ_clean,color="tab:cyan")
   l_fill_5        = fill_between(collect(4wsXPS.Nbe+1:5wsXPS.Nbe),H[4wsXPS.Nbe+1:5wsXPS.Nbe,:]*ρ_clean-σ_5*ones(Cdouble,wsXPS.Nbe),H[4wsXPS.Nbe+1:5wsXPS.Nbe,:]*ρ_clean+σ_5*ones(Cdouble,wsXPS.Nbe),alpha=0.5,color="tab:cyan")

   xlabel("SVD left index",fontsize=12)
   ylabel("PE signal [a.u.]",fontsize=12)
   legend([l_scat_1,(l_plot_1,l_fill_1),l_scat_2,(l_plot_2,l_fill_2),l_scat_3,(l_plot_3,l_fill_3),l_scat_4,(l_plot_4,l_fill_4),l_scat_5,(l_plot_5,l_fill_5)],["data 1"; "reconstruction 1"; "data 2"; "reconstruction 2"; "data 3"; "reconstruction 3"; "data 4"; "reconstruction 4"; "data 5"; "reconstruction 5"])






   ##
   ## cp
   ##



   A = [H; Matrix{Cdouble}(I,N,N); D_2nd; B];
   W_stop = ones(Cdouble,N);
   τ0 = 1.0e1 # 4
   x0 = zeros(Cdouble,N);
   s0 = A*x0;
   N_max_iter = 2*1000#0 # 00;
   r_n_tol=1.0
   r_y_tol=0.005;
   γ_H = 1.0e-2H;

   # xn,sn,taun,N_last = alg2_cp_gaussian_un_no_mem(x0,s0,IA_1,A,γ,γ_H,γ_D,W_stop;tau0=τ0,Niter=N_max_iter,r_n_tol=r_n_tol,r_y_tol=r_y_tol);
   # xn,sn,taun,F_ALL,G_ALL,F_STAR_ALL,G_STAR_ALL,INNER_PROD,T_ALL,N_last= alg2_cp_gaussian_un_no_mem(x0,s0,IA_1,A,γ,γ_H,γ_D,W_stop;tau0=τ0,Niter=N_max_iter,r_n_tol=r_n_tol,r_y_tol=r_y_tol);

   Ns = 5 # 00;
   ρ_samples = zeros(Cdouble,Ns,N);

   ρ_samples[1,:],_,_,_,_,_,_,_,_,N_last= alg2_cp_gaussian_un_no_mem_val(x0,s0,H*ρA_1,ρ0,ρB,A,γ,γ_H,γ_D,γ0,γB,W_stop;tau0=τ0,Niter=N_max_iter,r_n_tol=r_n_tol,r_y_tol=r_y_tol);
   t_elapsed = @elapsed for i in 2:Ns
      idx = shuffle_data(wsXPS.Nke*wsXPS.Nbe,wsXPS.Nke;Nmin=Nmin);
      # A = [H[idx,:]; Matrix{Cdouble}(I,N,N); D_2nd; B];
      x0 = zeros(Cdouble,N);
      s0 = A*x0
      # ρ_samples[i,:],_,_,_,_,_,_,_,_,N_last= alg2_cp_gaussian_un_no_mem_val(x0,s0,IA_1[idx],ρ0,ρB,A,γ[idx],γ_H[idx,:],γ_D,γ0,γB,W_stop;tau0=τ0,Niter=N_max_iter,r_n_tol=r_n_tol,r_y_tol=r_y_tol);
      ρ_samples[i,:],_,_,_,_,_,_,_,_,N_last= alg2_cp_gaussian_un_no_mem_val(x0,s0,H*ρA_1+sqrt.(γ).*randn(wsXPS.Nke*wsXPS.Nbe),ρ0,ρB,A,γ,γ_H,γ_D,γ0,γB,W_stop;tau0=τ0,Niter=N_max_iter,r_n_tol=r_n_tol,r_y_tol=r_y_tol);
   end

   μ_ρ = dropdims(mean(ρ_samples,dims=1),dims=1);
   σ_ρ = dropdims(std(ρ_samples,dims=1),dims=1);






   fig1 = figure()
   # state reconstruction
   ax1 = subplot(121)
   l_scat_1               = scatter(Zi,ρA_1)
   # l_plot_1_gt,           = plot(Zi,F.V[:,1:rankK]*F.Vt[1:rankK,:]*ρA_1[idx_res],color="tab:blue");
   l_plot_1_nso_un,       = plot(Zi,ρ_samples[1,:],color="tab:orange")
   l_plot_1_nso_un_sam,   = plot(Zi,μ_ρ,color="tab:pink")
   l_plot_1_nso_un_sam_σ  = fill_between(Zi,μ_ρ-σ_ρ,μ_ρ+σ_ρ,alpha=0.5,color="tab:pink")
   # l_plot_1_nso_un_sam_σ  = fill_between(Zi,μ_ρ-2.0σ_ρ,μ_ρ+2.0σ_ρ,alpha=0.5,color="tab:pink")
   xlabel("depth [nm]",fontsize=12)
   ylabel("concentration [a.u.]",fontsize=12)
   # legend([l_scat_1,l_plot_1_gt,l_plot_1_nso_un,(l_plot_1_nso_un_sam,l_plot_1_nso_un_sam_σ)],["GT","GT low rank","NSO model uncertainty","NSO sampling"])
   legend([l_scat_1,l_plot_1_nso_un,(l_plot_1_nso_un_sam,l_plot_1_nso_un_sam_σ)],["GT","CP model uncertainty","CP sampling"])


   # reconstructed data
   ax2 = subplot(122)
   # l_scat_2               = scatter(Ke,IA_1_clean,color="tab:green")
   l_scat_2_data          = scatter(wsXPS.Ke,H*ρA_1,color="tab:blue")
   l_plot_2_nso_un,       = plot(wsXPS.Ke,H*ρ_samples[1,:],color="tab:orange")
   l_plot_2_nso_un_sam,   = plot(wsXPS.Ke,H*μ_ρ,color="tab:pink")
   l_plot_2_nso_un_sam_σ  = fill_between(wsXPS.Ke,H*μ_ρ-sqrt.(γ),H*μ_ρ+sqrt.(γ),alpha=0.5,color="tab:pink")
   xlabel("kinetic energy [eV]",fontsize=12)
   ylabel("PE signal [a.u.]",fontsize=12)
   # legend([l_scat_2_data,l_scat_2,l_plot_2_nso_un,(l_plot_2_nso_un_sam,l_plot_2_nso_un_sam_σ)],["data","GT","NSO model uncertainty","NSO sampling"])
   legend([l_scat_2_data,l_plot_2_nso_un,(l_plot_2_nso_un_sam,l_plot_2_nso_un_sam_σ)],["data","CP model uncertainty","CP sampling"])


   fig1.set_figwidth(10.0) # (5.6)
   fig1.set_figheight(4.5)
   # tight_layout(pad=1.0, w_pad=0.2, h_pad=0.2)
   tight_layout(pad=1.0, w_pad=0.5, h_pad=0.5)

   ax1.annotate("a)", xy=(3, 1),  xycoords="data", xytext=(-0.1, 0.97), textcoords="axes fraction", color="black",fontsize=14)
   ax1.annotate("b)", xy=(3, 1),  xycoords="data", xytext=(1.02, 0.97), textcoords="axes fraction", color="black",fontsize=14)
end
