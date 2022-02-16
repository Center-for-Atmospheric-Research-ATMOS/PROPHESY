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
pebl_dict = Dict{Int64,Array{Cdouble,2}}();
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
   # plot(data[:,3],data[:,1]-data[:,2])
   plot(data[:,3],data[:,1])
   plot(data[:,3],data[:,2])
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



#TODO: plot normalied area w.r.t. depth

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

ANC1s[1:4]./ANO1s[1:4]

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
H_dict     = Dict{Int64,Tuple{Array{Cdouble,2},Array{Cdouble,2}}}();
for (i,data) in CSs
   Ke = pebl_dict[i][:,3]
   # x_dist_1 = data[4]*(τt[i,1]/(sqrt(2pi)*σt[i,1]))*exp.(-0.5*((Ke.-μt[i,1])/σt[i,1]).^2);
   # x_dist_2 = data[4]*(τt[i,2]/(sqrt(2pi)*σt[i,2]))*exp.(-0.5*((Ke.-μt[i,2])/σt[i,2]).^2);
   # # println(sum(x_dist_1+x_dist_2)*data_dict[i].Eb_step_eV[1])
   # figure(i)
   # plot(Ke,data[4]*data[1][:,1])
   # plot(Ke,(pebl_dict[i][:,1]-pebl_dict[i][:,2]))
   # fill_between(Ke,data[4]*(data[1][:,1]-sqrt.(diag(data[3]))),data[4]*(data[1][:,1]+sqrt.(diag(data[3]))),alpha=0.5,color="tab:blue")
   # plot(Ke,x_dist_1)
   # plot(Ke,x_dist_2)
   # plot(Ke,x_dist_1+x_dist_2)

   # cross section
   σ1s_peak1 = zeros(Cdouble,1,length(Ke));
   σ1s_peak2 = zeros(Cdouble,1,length(Ke));
   σ1s_peak1[1,:] = data_dict[i].cross_section[1]*(τt[i,1]/(sqrt(2pi)*σt[i,1]))*exp.(-0.5*((Ke.-μt[i,1])/σt[i,1]).^2);
   σ1s_peak2[1,:] = data_dict[i].cross_section[1]*(τt[i,2]/(sqrt(2pi)*σt[i,2]))*exp.(-0.5*((Ke.-μt[i,2])/σt[i,2]).^2);
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
   Be_exp[1,:] = ħν_exp.-Ke;
   setindex!(XPS_peak,(XPSsetup(ħν_exp,Fν_exp,μKe_exp,T_exp,Be_exp,σ1s_peak1;α_exp=1.0),XPSsetup(ħν_exp,Fν_exp,μKe_exp,T_exp,Be_exp,σ1s_peak2;α_exp=1.0)),i) #WARNING: not the right depths!!!!!!
   # wsXPS_1 = XPSsetup(ħν_exp,Fν_exp,μKe_exp,T_exp,Be_exp,σ1s_peak1;α_exp=1.0);
   # create the measurement model
   setindex!(H_dict,(Ψ_lin_peaks(Zi,XPS_peak[i][1];Nz=Nz0,σ_z=σ_z0,z0=z00,κ_cs=0.0,κ_eal=0.0),Ψ_lin_peaks(Zi,XPS_peak[i][2];Nz=Nz0,σ_z=σ_z0,z0=z00,κ_cs=0.0,κ_eal=0.0)),i);
   figure(i)
   subplot(121); imshow(H_dict[i][1]); colorbar(); subplot(122); imshow(H_dict[i][2]); colorbar()
end
