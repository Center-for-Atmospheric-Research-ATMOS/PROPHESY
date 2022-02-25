## load the packages used in the estimation
# plotting
using PyPlot
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
using XPSinv


##
## load some data: C1s and O1s
##

data_folder = "./data/HIPPIE_SDS_5mM/"

name_dict = Dict{Int64,String}();

setindex!(name_dict,string(data_folder,"SDS_5m_C1s_650eV.csv"),1)
setindex!(name_dict,string(data_folder,"SDS_5m_C1s_907eV.csv"),2)
setindex!(name_dict,string(data_folder,"SDS_5m_C1s_1315eV.csv"),3)
setindex!(name_dict,string(data_folder,"SDS_5m_C1s_1884eV.csv"),4)
setindex!(name_dict,string(data_folder,"SDS_5m_C1s_2200eV.csv"),5)

setindex!(name_dict,string(data_folder,"SDS_5m_O1s_897eV.csv"),6)
setindex!(name_dict,string(data_folder,"SDS_5m_O1s_1154eV.csv"),7)
setindex!(name_dict,string(data_folder,"SDS_5m_O1s_1565eV.csv"),8)
setindex!(name_dict,string(data_folder,"SDS_5m_O1s_2091eV.csv"),9)

data_dict = Dict{Int64,DataFrame}();

for (i,filename) in name_dict
   # load dataframe
   setindex!(data_dict,CSV.File(filename; delim=",", header=true) |> DataFrame,i)
   println(filename)
end

##
## baseline removal
##
pebl_dict = Dict{Int64,Array{Cdouble,2}}(); # [PE; BL; Eb]
for (i,data) in data_dict
   Ny = length(data[:,1]);
   D_2nd = diagm(Ny-2,Ny,1 => 2ones(Cdouble,Ny-2), 0 => -ones(Cdouble,Ny-2) ,2 => -ones(Cdouble,Ny-2));
   κb = 0.01;
   λb = 1.0e6 # 1.0e5;
   if i==1
      λb = 1.0e7;
   end
   z_baseline = baseline_removal(data[:,1],λb,κb,D_2nd);
   setindex!(pebl_dict,[reverse(data[:,1]) reverse(z_baseline) collect(data.Eb_start_eV[1]:data.Eb_step_eV[1]:data.Eb_end_eV[1])],i)
end


for (_,data) in pebl_dict
   figure();
   plot(data[:,3],data[:,1]-data[:,2])
   # plot(data[:,3],data[:,1])
   # plot(data[:,3],data[:,2])
end


# one Gaussian per peak: does not seem to fit well enough... or maybe it does... the residus are very similar in both cases
τmC = [0.9; 0.1]; # [0.75; 0.25]
μmC = [286.0; 288.0]; # [285.0; 287.0];
σmC = [0.5; 0.7]; # [285.0; 287.0]/500.0;

τmO = [1.0/2.0; 1.0/2.0];
μmO = [534.0; 536.0];
σmO = [0.7; 0.4] # [534.0; 536.0]/1000.0;

τt = zeros(Cdouble,9,2);
μt = zeros(Cdouble,9,2);
σt = zeros(Cdouble,9,2);
PLOT_FIG = true
for (i,data) in pebl_dict
   if i<=5
      # estimate the peaks centers and spreads
      if i==4
         idx_low = findfirst(data[:,3].>282.5);
      else
         idx_low = findfirst(data[:,3].>284.0);
      end
      if i==4
         idx_up  = findlast(data[:,3].<288.0);
      else
         idx_up  = findlast(data[:,3].<290.0);
      end
      # println(idx_low," ",idx_up," ",data[idx_low,3]," ",data[idx_up,3])
      # fit
      τt[i,:],μt[i,:],σt[i,:] = EM_peaks(data[idx_low:idx_up,3],data[idx_low:idx_up,1]-data[idx_low:idx_up,2],τmC,μmC,σmC,200)
   else
      idx_low = findfirst(data[:,3].>532.0);
      idx_up  = findlast(data[:,3].<538.0);
      Ke_shift = 0.0;
      if i==9
         idx_low = findfirst(data[:,3].>528.5);
         idx_up  = findlast(data[:,3].<535.0);
         Ke_shift = -2.0
      end
      if i==8
         idx_low = findfirst(data[:,3].>531.5);
         idx_up  = findlast(data[:,3].<537.5);
      end
      # fit
      if i==8
         τt[i,:],μt[i,:],σt[i,:] = EM_peaks(data[idx_low:idx_up,3],data[idx_low:idx_up,1]-data[idx_low:idx_up,2],τmO,μmO.+Ke_shift,σmO,200)
      else
         τt[i,:],μt[i,:],σt[i,:] = EM_peaks(data[idx_low:idx_up,3],data[idx_low:idx_up,1]-data[idx_low:idx_up,2],τmO,μmO.+Ke_shift,σmO,200)
      end
   end

   # component
   x_dist_1 = (τt[i,1]/σt[i,1])*exp.(-0.5*((data[:,3].-μt[i,1])/σt[i,1]).^2);
   x_dist_2 = (τt[i,2]/σt[i,2])*exp.(-0.5*((data[:,3].-μt[i,2])/σt[i,2]).^2);
   # plot
   if PLOT_FIG
      figure(i);
      plot(data[:,3],x_dist_1)
      plot(data[:,3],x_dist_2)
      plot(data[:,3],x_dist_1+x_dist_2)
      max_xx = maximum(x_dist_1+x_dist_2)
      max_dd = maximum(data[:,1]-data[:,2]);
      plot(data[:,3],(max_xx/max_dd)*(data[:,1]-data[:,2]))
   end
end




##
## cross section spread functions
##
σ_I1 = 0.1
σ_I2 = 0.05
σ_I3 = 0.1
σ_I4 = 0.05
σ_I5 = 0.01
σ_I6 = 0.2
σ_I7 = 0.1
σ_I8 = 0.05
σ_I9 = 0.05
σ_I = [σ_I1;σ_I2;σ_I3;σ_I4;σ_I5;σ_I6;σ_I7;σ_I8;σ_I9];

# number of iteration (main loop and inner loop of BFGS implementation)
Nbfgs = 10;
Nsearch = 10;
Nloop  = 10;

# estimate the cross section spread function and the overall integral for each spectrum)
CSs = Dict{Int64,Tuple{Array{Cdouble,2},Array{Cdouble,1},Array{Cdouble,2},Cdouble,Cdouble,Array{Cdouble,1}}}();
for (i,data) in pebl_dict
   if i==1
      setindex!(CSs,cross_section_spread_function_sample(data[:,1]-data[:,2],data[:,3],σ_I[i];Nbfgs=Nbfgs,Nsearch=Nsearch,N_sample=Nloop,σ_2nd=0.01),i)
   else
      setindex!(CSs,cross_section_spread_function_sample(data[:,1]-data[:,2],data[:,3],σ_I[i];Nbfgs=Nbfgs,Nsearch=Nsearch,N_sample=Nloop,σ_2nd=0.01),i)
   end
end

for (i,data) in CSs
   figure(i)
   plot(pebl_dict[i][:,3],data[1][:,1])
   plot(pebl_dict[i][:,3],(pebl_dict[i][:,1]-pebl_dict[i][:,2])/data[4])
   fill_between(pebl_dict[i][:,3],data[1][:,1]-sqrt.(diag(data[3])),data[1][:,1]+sqrt.(diag(data[3])),alpha=0.5,color="tab:blue")
end



κC = [7.31687612208259; 11.375; 10.2222222222222; 10.9368847079414]
τtC = [κC./(1.0.+κC) 1.0 ./(1.0.+κC)];
μtC = [285.7802 287.4658; 286.2876 288.0761; 286.1326 287.9222; 284.1892 285.9622]
σtC = [1.25360 1.25360; 1.21240 1.21240; 1.24480 1.24480; 1.32390 1.32390]

μtO = [534.4743 536.3316; 534.867 536.6754; 533.9868 535.7662; 531.3088 532.9999]
σtO = [1.6126 0.7004; 1.5767 0.8813; 1.6239 0.9968; 1.5338 1.1851]
areaO = [1448.862 5395.635; 4167.764 7380.07; 7125.726 5639.022; 1731.762 764.177]
κO = areaO[:,1]./areaO[:,2];
τtO = [κO./(1.0.+κO) 1.0 ./(1.0.+κO)]



if PLOT_FIG
   for (i,data) in pebl_dict
      if i<=4
         # component
         x_dist_1 = (τt[i,1]/σt[i,1])*exp.(-0.5*((data[:,3].-μt[i,1])/σt[i,1]).^2);
         x_dist_2 = (τt[i,2]/σt[i,2])*exp.(-0.5*((data[:,3].-μt[i,2])/σt[i,2]).^2);
         # other
         x_distC_1 = 2.0(τtC[i,1]/σtC[i,1])*exp.(-0.5*((data[:,3].-μtC[i,1])/0.5σtC[i,1]).^2);
         x_distC_2 = 2.0(τtC[i,2]/σtC[i,2])*exp.(-0.5*((data[:,3].-μtC[i,2])/0.5σtC[i,2]).^2);
         # plot
         figure(i);
         plot(data[:,3],x_dist_1,color="tab:blue")
         plot(data[:,3],x_dist_2,color="tab:blue")
         plot(data[:,3],x_dist_1+x_dist_2,color="tab:blue")
         plot(data[:,3],x_distC_1,color="tab:orange")
         plot(data[:,3],x_distC_2,color="tab:orange")
         plot(data[:,3],x_distC_1+x_distC_2,color="tab:orange")
         max_xx = maximum(x_dist_1+x_dist_2)
         max_dd = maximum(data[:,1]-data[:,2]);
         plot(data[:,3],(max_xx/max_dd)*(data[:,1]-data[:,2]),color="tab:green")
         max_xx = maximum(x_distC_1+x_distC_2)
         max_dd = maximum(data[:,1]-data[:,2]);
         plot(data[:,3],(max_xx/max_dd)*(data[:,1]-data[:,2]),color="tab:red")
      else
         if i>=6
            # component
            x_dist_1 = (τt[i,1]/σt[i,1])*exp.(-0.5*((data[:,3].-μt[i,1])/σt[i,1]).^2);
            x_dist_2 = (τt[i,2]/σt[i,2])*exp.(-0.5*((data[:,3].-μt[i,2])/σt[i,2]).^2);
            # other
            x_distO_1 = 2.0(τtO[i-5,1]/σtO[i-5,1])*exp.(-0.5*((data[:,3].-μtO[i-5,1])/0.5σtO[i-5,1]).^2);
            x_distO_2 = 2.0(τtO[i-5,2]/σtO[i-5,2])*exp.(-0.5*((data[:,3].-μtO[i-5,2])/0.5σtO[i-5,2]).^2);
            # plot
            figure(i);
            plot(data[:,3],x_dist_1,color="tab:blue")
            plot(data[:,3],x_dist_2,color="tab:blue")
            plot(data[:,3],x_dist_1+x_dist_2,color="tab:blue")
            plot(data[:,3],x_distO_1,color="tab:orange")
            plot(data[:,3],x_distO_2,color="tab:orange")
            plot(data[:,3],x_distO_1+x_distO_2,color="tab:orange")
            max_xx = maximum(x_dist_1+x_dist_2)
            max_dd = maximum(data[:,1]-data[:,2]);
            plot(data[:,3],(max_xx/max_dd)*(data[:,1]-data[:,2]),color="tab:green")
            max_xx = maximum(x_distO_1+x_distO_2)
            max_dd = maximum(data[:,1]-data[:,2]);
            plot(data[:,3],(max_xx/max_dd)*(data[:,1]-data[:,2]),color="tab:red")
         end
      end
   end
end


depthC1s = [data_dict[1].depth_nm[1]; data_dict[2].depth_nm[1]; data_dict[3].depth_nm[1]; data_dict[4].depth_nm[1]; data_dict[5].depth_nm[1]]
depthO1s = [data_dict[6].depth_nm[1]; data_dict[7].depth_nm[1]; data_dict[8].depth_nm[1]; data_dict[9].depth_nm[1]]

# FluxC1s = [data_dict[1].F[1]; data_dict[2].F[1]; data_dict[3].F[1]; -1.0; data_dict[5].F[1]]
FluxC1s = [data_dict[1].F[1]; data_dict[2].F[1]; data_dict[3].F[1]; data_dict[4].F[1]; data_dict[5].F[1]]
FluxO1s = [data_dict[6].F[1]; data_dict[7].F[1]; data_dict[8].F[1]; data_dict[9].F[1]]

CSC1s = [data_dict[1].cross_section[1]; data_dict[2].cross_section[1]; data_dict[3].cross_section[1]; data_dict[4].cross_section[1]; data_dict[5].cross_section[1]]
CSO1s = [data_dict[6].cross_section[1]; data_dict[7].cross_section[1]; data_dict[8].cross_section[1]; data_dict[9].cross_section[1]]

Ke_C1s = [data_dict[1].Eph_eV[1]-data_dict[1].Eb_eV[1]; data_dict[2].Eph_eV[1]-data_dict[2].Eb_eV[1]; data_dict[3].Eph_eV[1]-data_dict[3].Eb_eV[1]; data_dict[4].Eph_eV[1]-data_dict[4].Eb_eV[1]; data_dict[5].Eph_eV[1]-data_dict[5].Eb_eV[1]]
Ke_O1s = [data_dict[6].Eph_eV[1]-data_dict[6].Eb_eV[1]; data_dict[7].Eph_eV[1]-data_dict[7].Eb_eV[1]; data_dict[8].Eph_eV[1]-data_dict[8].Eb_eV[1]; data_dict[9].Eph_eV[1]-data_dict[9].Eb_eV[1]]

Eph_C1s = [data_dict[1].Eph_eV[1]; data_dict[2].Eph_eV[1]; data_dict[3].Eph_eV[1]; data_dict[4].Eph_eV[1]; data_dict[5].Eph_eV[1]]

AC1s = [CSs[1][4]; CSs[2][4]; CSs[3][4]; CSs[4][4]; CSs[5][4]]
AO1s = [CSs[6][4]; CSs[7][4]; CSs[8][4]; CSs[9][4]]

ANC1s = AC1s./(FluxC1s.*CSC1s)
ANO1s = AO1s./(FluxO1s.*CSO1s)


figure(); scatter(depthC1s[ANC1s.>0.0],ANC1s[ANC1s.>0.0]); ylim(0.0,1.1maximum(ANC1s[ANC1s.>0.0])) # 0.9minimum(ANC1s[ANC1s.>0.0])
figure(); scatter(depthO1s,ANO1s); ylim(0.0,1.1maximum(ANO1s))
# depthC1s[1:end-1]-depthO1s

# ANC1s[1:4]./ANO1s[1:4]

figure()
ledgend_id = []
for (i,data) in data_dict
   if i<=4
      global ledgend_id = [ledgend_id; i]
      plot(pebl_dict[i][:,3],(pebl_dict[i][:,1]-pebl_dict[i][:,2])/(data.F[1]*data.cross_section[1]))
   end
end



##
## create measurement model
##

## number of integration points for the Simpsons rule
Nz0 = 50;
## ρ_tot_int: total concentration. It should be mainly just water concentration
σ_z0 = 0.5; # [5 Å] width of the transition region
z00 = 0.0 # -0.5;  # half height depth

N = 50;
Z_max = 10.0;
Zi = collect(range(0.0,Z_max,length=N));
# Zi = collect(range(-1.0,Z_max,length=N)); #WARNING: in the model, the integral starts for 0, a negative depth would amplify the signal, which does not make sense, but would contribute to the signal! One may need to define the distance travel trhough the sample from far enough from the surface or set a saturation for negative depths


# CSs: Xend1,μ_XR1,Γ_XR1,R1,σ_R1,R_samples1
XPS_peak = Dict{Int64,Tuple{XPSsetup,XPSsetup}}();
XPS_peak1 = Dict{Int64,XPSsetup}();
XPS_peak2 = Dict{Int64,XPSsetup}();
H_dict     = Dict{Int64,Tuple{Array{Cdouble,2},Array{Cdouble,2}}}();
H_std_dict = Dict{Int64,Tuple{Array{Cdouble,2},Array{Cdouble,2}}}();
for (i,data) in CSs
   # emitted electron's kinetic energy
   Ke = data_dict[i].Eph_eV[1].-pebl_dict[i][:,3]

   # cross section
   σ1s_peak1 = zeros(Cdouble,1,length(Ke));
   σ1s_peak2 = zeros(Cdouble,1,length(Ke));
   σ1s_peak1[1,:] = data_dict[i].cross_section[1]*(τt[i,1]/(sqrt(2pi)*σt[i,1]))*exp.(-0.5*((pebl_dict[i][:,3].-μt[i,1])/σt[i,1]).^2);
   σ1s_peak2[1,:] = data_dict[i].cross_section[1]*(τt[i,2]/(sqrt(2pi)*σt[i,2]))*exp.(-0.5*((pebl_dict[i][:,3].-μt[i,2])/σt[i,2]).^2);
   # photon flux
   Fν_exp = [data_dict[i].F[1]]
   # photon energy
   ħν_exp = [1.0data_dict[i].Eph_eV[1]]
   # unknown transmission
   T_exp = [1.0];
   # central kinetic energy
   μKe_exp = [1.0(data_dict[i].Eph_eV[1]-data_dict[i].Eb_eV[1])];
   # binding energy
   Be_exp = zeros(Cdouble,1,length(Ke));
   Be_exp[1,:] = pebl_dict[i][:,3] # ħν_exp.-Ke;
   setindex!(XPS_peak,(XPSsetup(ħν_exp,Fν_exp,μKe_exp,T_exp,Be_exp,σ1s_peak1;α_exp=1.0),XPSsetup(ħν_exp,Fν_exp,μKe_exp,T_exp,Be_exp,σ1s_peak2;α_exp=1.0)),i) #WARNING: not the right depths!!!!!!
   setindex!(XPS_peak1,XPSsetup(ħν_exp,Fν_exp,μKe_exp,T_exp,Be_exp,σ1s_peak1;α_exp=1.0),i)
   setindex!(XPS_peak2,XPSsetup(ħν_exp,Fν_exp,μKe_exp,T_exp,Be_exp,σ1s_peak2;α_exp=1.0),i)
   XPS_peak[i][1].λe .= data_dict[i].depth_nm[1]
   XPS_peak[i][2].λe .= data_dict[i].depth_nm[1]
   XPS_peak1[i].λe   .= data_dict[i].depth_nm[1]
   XPS_peak2[i].λe   .= data_dict[i].depth_nm[1]
end


##
## peak area model
##

Hpeak1,Apeak1 = Ψ_lin_peaks_area(Zi,XPS_peak1;σ_z=σ_z0,z0=z00,κ_cs=0.05,κ_eal=0.05)
Hpeak2,Apeak2 = Ψ_lin_peaks_area(Zi,XPS_peak2;σ_z=σ_z0,z0=z00,κ_cs=0.05,κ_eal=0.05)

# normalize by the flux and the cross section
Apeak1_C1s = zeros(Cdouble,4,N);
Apeak1_C1s_std = zeros(Cdouble,4,N);
[Apeak1_C1s[i,:]     = Apeak1[i][1]/(data_dict[i].F[1]*data_dict[i].cross_section[1]) for i in 1:4];
[Apeak1_C1s_std[i,:] = Apeak1[i][2]/(data_dict[i].F[1]*data_dict[i].cross_section[1]) for i in 1:4];

Apeak2_C1s = zeros(Cdouble,4,N);
Apeak2_C1s_std = zeros(Cdouble,4,N);
[Apeak2_C1s[i,:]     = Apeak2[i][1]/(data_dict[i].F[1]*data_dict[i].cross_section[1]) for i in 1:4];
[Apeak2_C1s_std[i,:] = Apeak2[i][2]/(data_dict[i].F[1]*data_dict[i].cross_section[1]) for i in 1:4];

Apeak1_O1s = zeros(Cdouble,4,N);
Apeak1_O1s_std = zeros(Cdouble,4,N);
[Apeak1_O1s[i,:]     = Apeak1[i+5][1]/(data_dict[i+5].F[1]*data_dict[i+5].cross_section[1]) for i in 1:4];
[Apeak1_O1s_std[i,:] = Apeak1[i+5][2]/(data_dict[i+5].F[1]*data_dict[i+5].cross_section[1]) for i in 1:4];

Apeak2_O1s = zeros(Cdouble,4,N);
Apeak2_O1s_std = zeros(Cdouble,4,N);
[Apeak2_O1s[i,:]     = Apeak2[i+5][1]/(data_dict[i+5].F[1]*data_dict[i+5].cross_section[1]) for i in 1:4];
[Apeak2_O1s_std[i,:] = Apeak2[i+5][2]/(data_dict[i+5].F[1]*data_dict[i+5].cross_section[1]) for i in 1:4];



figure()
for i in 1:4
   plot(Zi,Apeak1_C1s[i,:])
   fill_between(Zi,Apeak1_C1s[i,:]-Apeak1_C1s_std[i,:],Apeak1_C1s[i,:]+Apeak1_C1s_std[i,:],alpha=0.5) #,color="tab:blue"
end

figure()
for i in 1:4
   plot(Zi,Apeak2_C1s[i,:])
   fill_between(Zi,Apeak2_C1s[i,:]-Apeak2_C1s_std[i,:],Apeak2_C1s[i,:]+Apeak2_C1s_std[i,:],alpha=0.5)
end


figure()
for i in 1:4
   plot(Zi,Apeak1_O1s[i,:])
   fill_between(Zi,Apeak1_O1s[i,:]-Apeak1_O1s_std[i,:],Apeak1_O1s[i,:]+Apeak1_O1s_std[i,:],alpha=0.5) #,color="tab:blue"
end

figure()
for i in 1:4
   plot(Zi,Apeak2_O1s[i,:])
   fill_between(Zi,Apeak2_O1s[i,:]-Apeak2_O1s_std[i,:],Apeak2_O1s[i,:]+Apeak2_O1s_std[i,:],alpha=0.5)
end

# check the structure of the null space: the basis that can be used to reconstruct the depth profile (that's what the data can say)
F_peak1_C1s = svd(Apeak1_C1s, full=true);
F_peak2_C1s = svd(Apeak2_C1s, full=true);

figure(); plot(Zi,F_peak1_C1s.Vt[1:4,:]')
figure(); plot(Zi,F_peak2_C1s.Vt[1:4,:]')


##
## relative peak area model
##

# check which peak to used for water (use the liquid water peak, not the gaz phase)
# now the data are the ratio between C1s peak area normalized by the photon flux and the cross section,
# and the O1s peak area normalized by the photon flux and the cross section
# Y = (A1C1s/Fσ)/(A1O1s/Fσ) : does not depend on α or T
y_all = ANC1s[1:4]./ANO1s;

# transform the peak area noise to normalized area noise
σ_all = zeros(Cdouble,9);
# [σ_all[i] = σ_I[i]*sqrt(size(pebl_dict[i],1))*data_dict[i].Eb_step_eV[1]/(data_dict[i].F[1]*data_dict[i].cross_section[1]) for i in 1:9]
[σ_all[i] = CSs[i][5]/(data_dict[i].F[1]*data_dict[i].cross_section[1]) for i in 1:9]

# relative area noise
σ_rel =  σ_all[1:4]./ANO1s # quite a good estimation of the noise level, I am happily surprised!
σ_rel[4] = 2.0σ_rel[4]; #

# that's the model! Use that to estimate the relative concentration (w.r.t. water concntration)
RApeak1 = Apeak1_C1s./depthO1s;         # it's also relative to the bulk concentration of water -> model for estimating ρA/ρwater
RApeak1_std = 0.01Apeak1_C1s_std./depthO1s;

D_1st = D1st(N);
D_2nd = D2nd(N);
B = zeros(Cdouble,2,N); B[1,1] = 1.0; B[2,end] = 1.0
γ = σ_rel.^2
Γ_inv = diagm(1.0./γ);
δ_D = 100000.0*(abs.(D_1st*logistic.(Zi.-z00,0.0,1.0,σ_z0)) + abs.(D_1st*logistic.(Zi.-1.0,0.0,1.0,σ_z0)) + abs.(D_1st*logistic.(Zi.-3.0,0.0,1.0,σ_z0)));
ι_D = 100000.0*(abs.(D_2nd*logistic.(Zi.-z00,0.0,1.0,σ_z0)) + abs.(D_2nd*logistic.(Zi.-1.0,0.0,1.0,σ_z0)) + abs.(D_2nd*logistic.(Zi.-3.0,0.0,1.0,σ_z0)));
ι_D .= 0.5;
ι_D .= 0.005;
γ_D = (ι_D .+ 0.01maximum(ι_D)).^2;
Γ_D_inv = diagm(1.0./γ_D);
ρ0 = 0.5;
γ0 = 1.0^2 # 0.01^2
ρB = 0.5;
γB = 1.0^2
ρ_bar = [ρ0;ρB];
γ_ρ = [γ0;γB];
Γ_ρ_inv = diagm(1.0./γ_ρ);
rankK = rank(RApeak1);


ρ_clean_nso,_,F = iterative_nso(rankK,RApeak1,Γ_inv,D_2nd,Γ_D_inv,B,Γ_ρ_inv,y_all,ρ_bar);

Ns = 50
ρ_samples_nso_h = zeros(Cdouble,Ns,N);
for j in 1:Ns
   ρ_samples_nso_h[j,:],_,_ = iterative_nso(rankK,RApeak1,Γ_inv,D_2nd,Γ_D_inv,B,Γ_ρ_inv,y_all+γ.*randn(4),ρ_bar);
end
μ_ρ_nso_h = dropdims(mean(ρ_samples_nso_h,dims=1),dims=1)
σ_ρ_nso_h = dropdims(std(ρ_samples_nso_h,dims=1),dims=1) # NOTE: the data are apparently good enough because the uncertainty in the estimation born by the data is fairly small, the model uncertainty is probably to be blamed for the bad reconstruction

# figure(); plot(Zi,ρ_clean_nso)


AA = [RApeak1; Matrix{Cdouble}(I,N,N); D_2nd; B];
W_stop = ones(Cdouble,N);
τ0 = 1.0e-10 # 1.0e-10
x0 = zeros(Cdouble,N) # ρA_1[idx_res]
s0 = AA*x0
N_max_iter = 200000 # 0;
r_n_tol=0.1*1.0
r_y_tol=0.1*0.005;


Ns = 50 # 0;
ρ_samples = zeros(Cdouble,Ns,N);

ρ_samples[1,:],_,_,_,_,_,_,_,_,N_last= alg2_cp_gaussian_un_no_mem_val(x0,s0,y_all,ρ0,ρB,AA,γ,RApeak1_std,γ_D,γ0,γB,W_stop;tau0=τ0,Niter=N_max_iter,r_n_tol=r_n_tol,r_y_tol=r_y_tol);
t_elapsed = @elapsed for i in 2:Ns
   # idx = shuffle_data(Nke,Npeak;Nmin=Nmin);
   # A = [H[idx,:]; Matrix{Cdouble}(I,N,N); D_2nd; B];
   x0 = zeros(Cdouble,N)
   s0 = AA*x0
   ρ_samples[i,:],_,_,_,_,_,_,_,_,N_last= alg2_cp_gaussian_un_no_mem_val(x0,s0,y_all+sqrt.(γ).*randn(4),ρ0,ρB,AA,γ,RApeak1_std,γ_D,γ0,γB,W_stop;tau0=τ0,Niter=N_max_iter,r_n_tol=r_n_tol,r_y_tol=r_y_tol);
end

μ_ρ = dropdims(mean(ρ_samples,dims=1),dims=1);
σ_ρ = dropdims(std(ρ_samples,dims=1),dims=1);


figure();
l_plot,        = plot(Zi,ρ_samples[1,:],color="tab:blue");
l_plot_μ,      = plot(Zi,μ_ρ,color="tab:pink");
l_fill_est     = fill_between(Zi,μ_ρ-σ_ρ,μ_ρ+σ_ρ,alpha=0.5,color="tab:pink");
xlabel("depth [nm]");
ylabel("normalized concentration \$\\frac{\\rho_{\\mathrm{SDS}}}{\\rho_{\\mathrm{bulk\\, water}}}\$");
xlim(Zi[1],Zi[end]);
ylim(-0.005)
legend([l_plot,(l_plot_μ,l_fill_est)],["CP estimation","data sampling: CP estimation"])
# savefig("SDS_concentration_profile_from_C1sO1s.png")
# savefig("SDS_concentration_profile_from_C1sO1s.pdf")

figure()
# scatter(y_all,RApeak1*ρ_clean_nso)
plot(y_all,y_all)
scatter(y_all,RApeak1*ρ_samples[1,:],color="tab:blue")
scatter(y_all,RApeak1*μ_ρ,color="tab:pink")
legend(["1:1","CP estimation","data sampling: CP estimation"])
xlim(0.015)
ylim(0.015)
xlabel("data")
ylabel("reconstructed data")
# savefig("SDS_data_vs_reconstruction_from_C1sO1s.png")
# savefig("SDS_data_vs_reconstruction_from_C1sO1s.pdf")


##
## simulated data
##

ρA_1 = logistic.(Zi.-2.0,0.0,1.0,2.0);
ρA_2 = logistic.(Zi.-2.0,0.0,1.0,2.0) .+ 2.0exp.(-(Zi.-1.0).^2. /(2.0*0.25^2));
ρA_3 = logistic.(Zi.-2.0,0.0,1.0,2.0) .+ exp.(-(Zi.-1.5).^2. /(2.0*0.5^2));
ρA_4 = exp.(-(Zi.-2.5).^2. /(2.0*0.5^2));


y_all_1 = RApeak1*ρA_1
y_all_2 = RApeak1*ρA_2
y_all_3 = RApeak1*ρA_3
y_all_4 = RApeak1*ρA_4



Ns = 50 # 0;
ρ_samples_1 = zeros(Cdouble,Ns,N);
ρ0_1 = 0.035; γ0_1 = 0.001^2
ρB_1 = 1.0; γB_1 = 0.001^2
γ_1 = (1.0e-4)^2*ones(Cdouble,4);
N_max_iter = 200000
γ_D_1 = (0.05(1.0e-6 .+D_2nd*ρA_1)).^2 # 0.0000001*γ_D
ρ_samples_1[1,:],_,_,_,_,_,_,_,_,N_last= alg2_cp_gaussian_un_no_mem_val(x0,s0,y_all_1,ρ0_1,ρB_1,AA,γ_1,RApeak1_std,γ_D_1,γ0_1,γB_1,W_stop;tau0=τ0,Niter=N_max_iter,r_n_tol=r_n_tol,r_y_tol=r_y_tol);
t_elapsed = @elapsed for i in 2:Ns
   # idx = shuffle_data(Nke,Npeak;Nmin=Nmin);
   # A = [H[idx,:]; Matrix{Cdouble}(I,N,N); D_2nd; B];
   x0 = zeros(Cdouble,N)
   s0 = AA*x0
   ρ_samples_1[i,:],_,_,_,_,_,_,_,_,N_last= alg2_cp_gaussian_un_no_mem_val(x0,s0,y_all_1+sqrt.(γ_1).*randn(4),ρ0_1,ρB_1,AA,γ_1,RApeak1_std,γ_D_1,γ0_1,γB_1,W_stop;tau0=τ0,Niter=N_max_iter,r_n_tol=r_n_tol,r_y_tol=r_y_tol);
end

μ_ρ_1 = dropdims(mean(ρ_samples_1,dims=1),dims=1);
σ_ρ_1 = dropdims(std(ρ_samples_1,dims=1),dims=1);


figure();
l_plot,        = plot(Zi,ρ_samples_1[1,:],color="tab:blue");
l_plot_μ,      = plot(Zi,μ_ρ_1,color="tab:pink");
l_fill_est     = fill_between(Zi,μ_ρ_1-σ_ρ_1,μ_ρ_1+σ_ρ_1,alpha=0.5,color="tab:pink");
l_plot_gt,     = plot(Zi,ρA_1,color="tab:green");
xlabel("depth [nm]");
ylabel("normalized concentration \$\\frac{\\rho_{\\mathrm{X}}}{\\rho_{\\mathrm{bulk\\, water}}}\$");
xlim(Zi[1],Zi[end]);
ylim(-0.005)
legend([l_plot,(l_plot_μ,l_fill_est),l_plot_gt],["CP estimation","data sampling: CP estimation","GT"])
# savefig("X1_concentration_profile_from_C1sO1s.png")
# savefig("X1_concentration_profile_from_C1sO1s.pdf")


figure()
plot(y_all_1,y_all_1)
scatter(y_all_1,RApeak1*ρ_samples_1[1,:],color="tab:blue")
legend(["1:1","CP estimation","data sampling: CP estimation"])
xlim(0.015)
ylim(0.015)
xlabel("data")
ylabel("reconstructed data")
# savefig("X1_data_vs_reconstruction_from_C1sO1s.png")
# savefig("X1_data_vs_reconstruction_from_C1sO1s.pdf")



Ns = 50 # 0;
ρ_samples_2 = zeros(Cdouble,Ns,N);
ρ0_2 = 0.035; γ0_2 = 0.001^2
ρB_2 = 1.0; γB_2 = 0.001^2
γ_2 = (1.0e-4)^2*ones(Cdouble,4);
N_max_iter = 200000
γ_D_2 = (0.05(1.0e-6 .+D_2nd*ρA_2)).^2
ρ_samples_2[1,:],_,_,_,_,_,_,_,_,N_last= alg2_cp_gaussian_un_no_mem_val(x0,s0,y_all_2,ρ0_2,ρB_2,AA,γ_2,RApeak1_std,γ_D_2,γ0_2,γB_2,W_stop;tau0=τ0,Niter=N_max_iter,r_n_tol=r_n_tol,r_y_tol=r_y_tol);
t_elapsed = @elapsed for i in 2:Ns
   # idx = shuffle_data(Nke,Npeak;Nmin=Nmin);
   # A = [H[idx,:]; Matrix{Cdouble}(I,N,N); D_2nd; B];
   x0 = zeros(Cdouble,N)
   s0 = AA*x0
   ρ_samples_2[i,:],_,_,_,_,_,_,_,_,N_last= alg2_cp_gaussian_un_no_mem_val(x0,s0,y_all_2+sqrt.(γ_2).*randn(4),ρ0_2,ρB_2,AA,γ_2,RApeak1_std,γ_D_2,γ0_2,γB_2,W_stop;tau0=τ0,Niter=N_max_iter,r_n_tol=r_n_tol,r_y_tol=r_y_tol);
end

μ_ρ_2 = dropdims(mean(ρ_samples_2,dims=1),dims=1);
σ_ρ_2 = dropdims(std(ρ_samples_2,dims=1),dims=1);


figure();
l_plot,        = plot(Zi,ρ_samples_2[1,:],color="tab:blue");
l_plot_μ,      = plot(Zi,μ_ρ_2,color="tab:pink");
l_fill_est     = fill_between(Zi,μ_ρ_2-σ_ρ_2,μ_ρ_2+σ_ρ_2,alpha=0.5,color="tab:pink");
l_plot_gt,     = plot(Zi,ρA_2,color="tab:green");
xlabel("depth [nm]");
ylabel("normalized concentration \$\\frac{\\rho_{\\mathrm{X}}}{\\rho_{\\mathrm{bulk\\, water}}}\$");
xlim(Zi[1],Zi[end]);
ylim(-0.005)
legend([l_plot,(l_plot_μ,l_fill_est),l_plot_gt],["CP estimation","data sampling: CP estimation","GT"])
# savefig("X2_concentration_profile_from_C1sO1s.png")
# savefig("X2_concentration_profile_from_C1sO1s.pdf")

figure()
plot(y_all_2,y_all_2)
scatter(y_all_2,RApeak1*ρ_samples_2[1,:],color="tab:blue")
legend(["1:1","CP estimation","data sampling: CP estimation"])
xlim(0.015)
ylim(0.015)
xlabel("data")
ylabel("reconstructed data")
# savefig("X2_data_vs_reconstruction_from_C1sO1s.png")
# savefig("X2_data_vs_reconstruction_from_C1sO1s.pdf")




Ns = 50 # 0;
ρ_samples_3 = zeros(Cdouble,Ns,N);
ρ0_3 = 0.035; γ0_3 = 0.001^2
ρB_3 = 1.0; γB_3 = 0.001^2
γ_3 = (1.0e-4)^2*ones(Cdouble,4); # (1.0e-3)^2*ones(Cdouble,4); #
N_max_iter = 200000
γ_D_3 = (0.05(1.0e-6 .+D_2nd*ρA_3)).^2
ρ_samples_3[1,:],_,_,_,_,_,_,_,_,N_last= alg2_cp_gaussian_un_no_mem_val(x0,s0,y_all_3,ρ0_3,ρB_3,AA,γ_3,RApeak1_std,γ_D_3,γ0_3,γB_3,W_stop;tau0=τ0,Niter=N_max_iter,r_n_tol=r_n_tol,r_y_tol=r_y_tol);
t_elapsed = @elapsed for i in 2:Ns
   # idx = shuffle_data(Nke,Npeak;Nmin=Nmin);
   # A = [H[idx,:]; Matrix{Cdouble}(I,N,N); D_2nd; B];
   x0 = zeros(Cdouble,N)
   s0 = AA*x0
   ρ_samples_3[i,:],_,_,_,_,_,_,_,_,N_last= alg2_cp_gaussian_un_no_mem_val(x0,s0,y_all_3+sqrt.(γ_3).*randn(4),ρ0_3,ρB_3,AA,γ_3,RApeak1_std,γ_D_3,γ0_3,γB_3,W_stop;tau0=τ0,Niter=N_max_iter,r_n_tol=r_n_tol,r_y_tol=r_y_tol);
end

μ_ρ_3 = dropdims(mean(ρ_samples_3,dims=1),dims=1);
σ_ρ_3 = dropdims(std(ρ_samples_3,dims=1),dims=1);


figure();
l_plot,        = plot(Zi,ρ_samples_3[1,:],color="tab:blue");
l_plot_μ,      = plot(Zi,μ_ρ_3,color="tab:pink");
l_fill_est     = fill_between(Zi,μ_ρ_3-σ_ρ_3,μ_ρ_3+σ_ρ_3,alpha=0.5,color="tab:pink");
l_plot_gt,     = plot(Zi,ρA_3,color="tab:green");
xlabel("depth [nm]");
ylabel("normalized concentration \$\\frac{\\rho_{\\mathrm{X}}}{\\rho_{\\mathrm{bulk\\, water}}}\$");
xlim(Zi[1],Zi[end]);
ylim(-0.005)
legend([l_plot,(l_plot_μ,l_fill_est),l_plot_gt],["CP estimation","data sampling: CP estimation","GT"])
# savefig("X3_concentration_profile_from_C1sO1s.png")
# savefig("X3_concentration_profile_from_C1sO1s.pdf")

figure()
plot(y_all_3,y_all_3)
scatter(y_all_3,RApeak1*ρ_samples_3[1,:],color="tab:blue")
legend(["1:1","CP estimation","data sampling: CP estimation"])
xlim(0.015)
ylim(0.015)
xlabel("data")
ylabel("reconstructed data")
# savefig("X3_data_vs_reconstruction_from_C1sO1s.png")
# savefig("X3_data_vs_reconstruction_from_C1sO1s.pdf")


Ns = 50 # 0;
ρ_samples_4 = zeros(Cdouble,Ns,N);
ρ0_4 = 0.035; γ0_4 = 0.001^2
ρB_4 = 0.0; γB_4 = 0.001^2
γ_4 = (1.0e-4)^2*ones(Cdouble,4);
N_max_iter = 200000
γ_D_4 = (0.05(1.0e-6 .+D_2nd*ρA_4)).^2
ρ_samples_4[1,:],_,_,_,_,_,_,_,_,N_last= alg2_cp_gaussian_un_no_mem_val(x0,s0,y_all_4,ρ0_4,ρB_4,AA,γ_4,RApeak1_std,γ_D_4,γ0_4,γB_4,W_stop;tau0=τ0,Niter=N_max_iter,r_n_tol=r_n_tol,r_y_tol=r_y_tol);
t_elapsed = @elapsed for i in 2:Ns
   # idx = shuffle_data(Nke,Npeak;Nmin=Nmin);
   # A = [H[idx,:]; Matrix{Cdouble}(I,N,N); D_2nd; B];
   x0 = zeros(Cdouble,N)
   s0 = AA*x0
   ρ_samples_4[i,:],_,_,_,_,_,_,_,_,N_last= alg2_cp_gaussian_un_no_mem_val(x0,s0,y_all_4+sqrt.(γ_4).*randn(4),ρ0_4,ρB_4,AA,γ_4,RApeak1_std,γ_D_4,γ0_4,γB_4,W_stop;tau0=τ0,Niter=N_max_iter,r_n_tol=r_n_tol,r_y_tol=r_y_tol);
end

μ_ρ_4 = dropdims(mean(ρ_samples_4,dims=1),dims=1);
σ_ρ_4 = dropdims(std(ρ_samples_4,dims=1),dims=1);

figure();
l_plot,        = plot(Zi,ρ_samples_4[1,:],color="tab:blue");
l_plot_μ,      = plot(Zi,μ_ρ_4,color="tab:pink");
l_fill_est     = fill_between(Zi,μ_ρ_4-σ_ρ_4,μ_ρ_4+σ_ρ_4,alpha=0.5,color="tab:pink");
l_plot_gt,     = plot(Zi,ρA_4,color="tab:green");
xlabel("depth [nm]");
ylabel("normalized concentration \$\\frac{\\rho_{\\mathrm{X}}}{\\rho_{\\mathrm{bulk\\, water}}}\$");
xlim(Zi[1],Zi[end]);
ylim(-0.005)
legend([l_plot,(l_plot_μ,l_fill_est),l_plot_gt],["CP estimation","data sampling: CP estimation","GT"])
# savefig("X4_concentration_profile_from_C1sO1s.png")
# savefig("X4_concentration_profile_from_C1sO1s.pdf")

figure()
plot(y_all_4,y_all_4)
scatter(y_all_4,RApeak1*ρ_samples_4[1,:],color="tab:blue")
legend(["1:1","CP estimation","data sampling: CP estimation"])
xlim(0.015)
ylim(0.015)
xlabel("data")
ylabel("reconstructed data")
# savefig("X4_data_vs_reconstruction_from_C1sO1s.png")
# savefig("X4_data_vs_reconstruction_from_C1sO1s.pdf")
