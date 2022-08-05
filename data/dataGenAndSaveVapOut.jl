## load the packages used in the estimation
# plotting
using PyPlot
rc("text", usetex=true)
rc("figure",max_open_warning=50)
using myPlot

# data manipulation (loading, writing, etc)
using Printf
using XLSX # CSV does not deal with multiple sheets
using DataFrames

# scientific package from the official Julia repositories
using LinearAlgebra
using StatsBase
using Interpolations

# implemented scientific packages
using utilsFun  # for the softMax functions

# modeling XPS
using XPSpack # experiment model (geometry factor and cross section estimation)
using ATTIRE  # kinetic energy analyzer



# tags
SHORT_RANGE = true              # select either wide range of attenuation lengths (false) or a restricted range more similar to experimental setup (true)

MODEL_5   = true               # select the number of attenuation lengths probed
MODEL_10  = false
MODEL_20  = false

N_model_sample = 100;           # 100 should be more than enough for 5 attenuation length, but maybe not for 20

FLAG_0001 = false               # selection of the profile (one must be true and the others false)
FLAG_0002 = true
FLAG_0003 = false
FLAG_0004 = false

save_folder = "./";
SAVE_DATA = false   # flag for saving data and model
SAVE_FIG  = false   # save the simulated data or not

# geometry setup
λe0 = 2.0e-3;        # reference penetration depth in μm
δr = 2.0e-3          # transition to vacuum layer thickness (let's set about 1 nm)
k0 = 10;             # compute the measurement model over a distance of k0*λe0
Nr = 201;            # number of discretization points in the radial dimension
Nr_lowres = 101 ;    # for the low resolution model
Nθ = 256;            # number of discretization points in the polar angle dimension
Ny = 256;            # number of discretization points in the cylinder axis dimension
μ0 = 10.0; # 20.0;   # radius of the cylinder
L = 25; # 100.0;     # height of the irradiated sample (the vertical extent of the beam is more like 20μm instead of 100μm)
x0 = sqrt(2.0)*500.0 # (x0,y0,z0) are the coordinates of the analyzer's apperture
y0 = 0.0;
z0 = 500.0;

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


if FLAG_0001
    ρA_1 = logistic.(1000.0*(μ0.-r)/0.5,ρ_vac,ρ0,1.0);
    exp_tag     = "0001"
end
if FLAG_0002
    ρA_1 = logistic.(1000.0*(μ0.-r)/0.5,ρ_vac,ρ0,1.0) .+ 2.0exp.(-(1000.0*(μ0.-r).-0.0).^2. /(2.0*0.25^2));
    exp_tag     = "0002"
end
if FLAG_0003
    ρA_1 = logistic.(1000.0*(μ0.-r)/0.5,ρ_vac,ρ0,1.0) .+ exp.(-(1000.0*(μ0.-r).-0.5).^2. /(2.0*0.5^2));
    exp_tag     = "0003"
end
if FLAG_0004
    ρA_1 = exp.(-(1000.0((μ0.-r))).^2. /(2.0*0.5^2));
    exp_tag     = "0004"
end

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
    λe1 = 1.3; hν1 = 365.0;
    λe2 = 2.5; hν2 = 1500.0;
    save_folder = string(save_folder,"eal_",Ndata,"_restricted_range/")
else
    λe1 = 0.5; hν1 = 310.0;
    λe2 = 5.5; hν2 = 1900.0;
    save_folder = string(save_folder,"eal_",Ndata,"/")
end
λe = 1.0e-3collect(range(λe1,λe2,Ndata));              # attenuation length range
hν = collect(LinRange(hν1, hν2,Ndata));                # central photon energy for each measurement

##
## alignment
##
xc = 100.0 # 200.0 # 400.0;
yc = 99.0 # 75.0; # 0.0;
σx = 100.0;
σy = 25.0;
bp = beamProfile(xc,yc,σx,σy);
α_al = zeros(Cdouble,Ndata);
for i in 1:Ndata
    H_r,H_rθy,H_r_ph,H_rθy_ph,_,_,_,α_al[i] =  alignmentParameter(bp,r,θ,y,x0,y0,z0,μ0,λe[i]);
end

## 
## cross section density model: just a dummy model looking like C1s
## 

dKe = 0.05;
Be = collect(286.0:dKe:298.0);
Nspectrum = length(Be);
μBe = [290.2; 292.0; 293.0]
σ_be = sqrt(2.0)*[0.45; 0.25; 0.6];
BeC1s = mean(Be);
ΔBeC1s = (Be[end]-Be[1])/5;


function σ_cs(hν::Cdouble,Ke::Cdouble,μKe::Cdouble;μKe0::Cdouble=50.0,μKe1::Cdouble=1200.0)
    Be = hν-Ke;
    # partial cross section (one for each chemical state)
    σ_peak_1 = (1.0/sqrt(2.0π*σ_be[1]^2))*exp.(-(Be-μBe[1])^2/(2.0σ_be[1]^2));
    σ_peak_2 = (1.0/sqrt(2.0π*σ_be[2]^2))*exp.(-(Be-μBe[2])^2/(2.0σ_be[2]^2));
    σ_peak_3 = (1.0/sqrt(2.0π*σ_be[3]^2))*exp.(-(Be-μBe[3])^2/(2.0σ_be[3]^2));
    # quantity of chemical states
    p1 = 0.85 .+ (0.77-0.85)*(μKe.-μKe0)./(μKe1-μKe0);
    p2 = 0.125 .+ (0.12-0.125)*(μKe.-μKe0)./(μKe1-μKe0);
    p3 = 1.0-(p1+p2);

    # cross section value (only for hν ∈ [295,1500.0])
    XPSpack.σ_C1s_interp[hν]*(p1*σ_peak_1+p2*σ_peak_2+p3*σ_peak_3)
end

##
## acqusition parameters
##

θ_aperture = 0.5*π/4                                                       # aperture's cone angle
α_Ω        = 4π*sin(θ_aperture/2.0)^2;                                     # aperture's solid angle
Tj         = α_Ω*(10.0.+0.0collect(LinRange(5.0,10.0,Ndata)));             # transmission factors
σ_ke       = 2.0*dKe*ones(Cdouble,Ndata);                                  # kinetic energy bandwidths of the analyzer (one per photon energy)
dhν        = hν.*((1.0/25000.0)*(hν.<500.0) + (1.0/15000.0)*(hν.>=500.0)); # bandwidth of the photon beam
Fνj        = 1.5e11*ones(Cdouble,Ndata);                                   # flux densities


# dictionary where to push the data and geometry factor
dictAllData = Dict()
dictAllGeom = Dict()
for j in 1:Ndata # can potentially multi-thread this loop: but need to sync before writing files
    # for a given photon energy νj, measure a spectrum
    local Keij = reverse(hν[j] .- Be) ;                                                        # centers of the analyzer's channels
    local kKe = floor(5σ_ke[j]/dKe);                                                           # number of extra discretization point so that the channels at Ki[1] and Ki[end] are not missing to much information
    local Kdisc = [Keij[1].-dKe*reverse(collect(1:kKe)); Keij; Keij[end].+dKe*collect(1:kKe)]; # discretization point
    local μKe = 0.5*(Keij[1]+Keij[end]);                                                       # central kinetic energy (a bit the same role as pass energy)
    
    ##
    ## compute the cross section of the sample (to be estimated from the data in an estimation setting)
    ##

    local Be0 = hν[j] .- Keij;
    local σ_cs_fg =  σ_cs.(hν[j],Keij,μKe);

    ##
    ## electron flux without the geometry factor: signal of interest (Sj) and background signal (Sj_bg)
    ##
    local Sj,Sj_bg,_,_ = simulateSpectrum(Fνj[j],hν[j],dhν[j],
        Keij,σ_ke[j],Tj[j],Kdisc,
        reverse(Be0),reverse(σ_cs_fg),σ_bg(μKe)*σ_bg_density(Keij,hν[j]-BeC1s,ΔBeC1s));

    ##
    ## geometry factors
    ##
    local H_deom,_,H_geom,_,_,_,_,_ = alignmentParameter(bp,r,θ,y,x0,y0,z0,μ0,λe[j])

    ##
    ## for the signal without noise
    ##
    local SbgC1s      = (H_geom'*ρA_1)*Sj_bg
    local SpectrumA_1 = (H_geom'*ρA_1)*Sj

    ##
    ## add noise
    ##
    local SC1snoise = countElectrons(SbgC1s+SpectrumA_1)

    # plot signals w.r.t. the kinetic energy
    # figure(); plot(Keij,SbgC1s); plot(Keij,SbgC1s+SpectrumA_1); scatter(Keij,rand.(Poisson.(SbgC1s+SpectrumA_1))); xlabel("kinetic energy [eV]"); ylabel("spectrum [a.u.]") 

    ##
    ## push data to dicts
    ##

    local dictData = Dict( "Ke" => Keij, "Be" => Be0, "μKe" => μKe,
        "σ_cs_dens" => σ_cs_fg./XPSpack.σ_C1s_interp[hν[j]], "σ_tot" => XPSpack.σ_C1s_interp[hν[j]], 
        "SpectrumA_1" => SpectrumA_1, "Sbg" => SbgC1s, 
        "Stot" => SbgC1s+SpectrumA_1, "Snoisy" => SC1snoise,
        "T" => Tj[j], "λ" => 1.0e3λe[j], "F" => Fνj[j], "hν" => hν[j]);

    local dictGeom = Dict("model" => "sharp edge cylinder + outside vapor", 
                "hν" => hν[j], "λ" => 1.0e3λe[j], "radius" => μ0, "max_depth" => k0*λe0,
                    "x0" => x0, "y0" => y0, "z0" => z0, "δr" => δr,
                    "r" => r, "H" => H_deom, "H_true" => H_geom, "ρ" => ρA_1)

    local dfGeom = DataFrame(dictGeom);
    local dfData = DataFrame(dictData);

    dictAllData[Symbol(string("hν_",Int64(round(hν[j]))))] = (eachcol(dfData),names(dfData))
    dictAllGeom[Symbol(string("λe_",string(1.0e3λe[j])))]  = (eachcol(dfGeom),names(dfGeom))
end

if SAVE_DATA
    mkpath(string(save_folder,exp_tag,"/"));
    XLSX.writetable(string(save_folder,exp_tag,"/data.xlsx"); dictAllData...) # TODO: get outside the loop
    XLSX.writetable(string(save_folder,exp_tag,"/model.xlsx"); dictAllGeom...)
end
include("plotData.jl")


# TODO: simulate data for several photon energy and fit the peaks in order to check how they shift, SpectrumA_1 w.r.t. reverse(Be)

# [a] E. Kukk, J.D. Bozek, G. Snell, W.-T. Cheng and N. Berrah, Phys. Rev. A v.63, 062702 (2001).
# [b] E. Kukk, K. Ueda, U. Hergenhahn, J. Liu X, G. Prumper, H. Yoshida, Y. Tamenori, C. Makochekanwa, T. Tanaka, M. Kitajima and H. Tanaka, Phys.Rev.Lett. v. 95, p. 133001 (2005).

τm = [0.85; 0.125; 1.0-0.85-0.125];  # [1.0/3.0; 1.0/3.0; 1.0/3.0];
μm = [290.2; 292.0; 293.0]; # [290.3; 291.9; 293.5];
figure()
for j in  1:Ndata # 1:5:Ndata
    local μKe0=50.0
    local μKe1=1200.0
    local symbol_h = Symbol(string("hν_",Int64(round(hν[j]))));
    local be = dictAllData[symbol_h][1][:Be];
    local μKe = dictAllData[symbol_h][1][:μKe];
    # partial cross section (one for each chemical state)
    σ_peak_1 = (1.0/sqrt(2.0π*σ_be[1]^2))*exp.(-(be.-μBe[1]).^2/(2.0σ_be[1]^2));
    σ_peak_2 = (1.0/sqrt(2.0π*σ_be[2]^2))*exp.(-(be.-μBe[2]).^2/(2.0σ_be[2]^2));
    σ_peak_3 = (1.0/sqrt(2.0π*σ_be[3]^2))*exp.(-(be.-μBe[3]).^2/(2.0σ_be[3]^2));
    # quantity of chemical states
    p1 = 0.85 .+ (0.77-0.85)*(μKe[1].-μKe0)./(μKe1-μKe0);
    p2 = 0.125 .+ (0.12-0.125)*(μKe[1].-μKe0)./(μKe1-μKe0);
    p3 = 1.0-(p1+p2);

    # estimation
    σ_est_1 = τt[j,1]*(1.0/sqrt(2.0π*σt[j,1]^2))*exp.(-(be.-μt[j,1]).^2/(2.0σt[j,1]^2));
    σ_est_2 = τt[j,2]*(1.0/sqrt(2.0π*σt[j,2]^2))*exp.(-(be.-μt[j,2]).^2/(2.0σt[j,2]^2));
    σ_est_3 = τt[j,3]*(1.0/sqrt(2.0π*σt[j,3]^2))*exp.(-(be.-μt[j,3]).^2/(2.0σt[j,3]^2));

    println(dKe*sum(p1*σ_peak_1+p2*σ_peak_2+p3*σ_peak_3))
    println(dKe*sum(σ_est_1+σ_est_2+σ_est_3),"\n")

    plot(be,p1*σ_peak_1+p2*σ_peak_2+p3*σ_peak_3)
    scatter(be,σ_est_1+σ_est_2+σ_est_3)

    plot(be,(dictAllData[symbol_h][1][:Snoisy]-dictAllData[symbol_h][1][:Sbg])/(dKe*sum(dictAllData[symbol_h][1][:Snoisy]-dictAllData[symbol_h][1][:Sbg])))
end

σm = sqrt(2.0)*[0.45; 0.25; 0.6]; # [290.3; 291.9; 293.5]/500.0;


τt = zeros(Cdouble,Ndata,3);
μt = zeros(Cdouble,Ndata,3);
σt = zeros(Cdouble,Ndata,3);
for j in 1:Ndata
    local symbol_h = Symbol(string("hν_",Int64(round(hν[j]))));
    local be = dictAllData[symbol_h][1][:Be]; # reverse();
    # local spectrum = dictAllData[symbol_h][1][:SpectrumA_1];
    local spectrum = dictAllData[symbol_h][1][:Snoisy]-dictAllData[symbol_h][1][:Sbg];
   # estimate the peaks centers and spreads
   τt[j,:],μt[j,:],σt[j,:] = EM_peaks(be,spectrum,τm,μm,σm,200)
end

figure(); 
scatter(collect(1:Ndata),abs.(μt[:,1].-μBe[1]))
scatter(collect(1:Ndata),abs.(μt[:,2].-μBe[2]))
scatter(collect(1:Ndata),abs.(μt[:,3].-μBe[3]))

figure(); 
scatter(collect(1:Ndata),σt[:,1])
scatter(collect(1:Ndata),σt[:,2])
scatter(collect(1:Ndata),σt[:,3])

figure(); 
scatter(collect(1:Ndata),τt[:,1])
scatter(collect(1:Ndata),τt[:,2])
scatter(collect(1:Ndata),τt[:,3])

figure()
for j in  1:Ndata # 1:5:Ndata
    local μKe0=50.0
    local μKe1=1200.0
    local symbol_h = Symbol(string("hν_",Int64(round(hν[j]))));
    local be = dictAllData[symbol_h][1][:Be];
    local μKe = dictAllData[symbol_h][1][:μKe];
    # partial cross section (one for each chemical state)
    σ_peak_1 = (1.0/sqrt(2.0π*σ_be[1]^2))*exp.(-(be.-μBe[1]).^2/(2.0σ_be[1]^2));
    σ_peak_2 = (1.0/sqrt(2.0π*σ_be[2]^2))*exp.(-(be.-μBe[2]).^2/(2.0σ_be[2]^2));
    σ_peak_3 = (1.0/sqrt(2.0π*σ_be[3]^2))*exp.(-(be.-μBe[3]).^2/(2.0σ_be[3]^2));
    # quantity of chemical states
    p1 = 0.85 .+ (0.77-0.85)*(μKe[1].-μKe0)./(μKe1-μKe0);
    p2 = 0.125 .+ (0.12-0.125)*(μKe[1].-μKe0)./(μKe1-μKe0);
    p3 = 1.0-(p1+p2);

    # estimation
    σ_est_1 = τt[j,1]*(1.0/sqrt(2.0π*σt[j,1]^2))*exp.(-(be.-μt[j,1]).^2/(2.0σt[j,1]^2));
    σ_est_2 = τt[j,2]*(1.0/sqrt(2.0π*σt[j,2]^2))*exp.(-(be.-μt[j,2]).^2/(2.0σt[j,2]^2));
    σ_est_3 = τt[j,3]*(1.0/sqrt(2.0π*σt[j,3]^2))*exp.(-(be.-μt[j,3]).^2/(2.0σt[j,3]^2));

    println(dKe*sum(p1*σ_peak_1+p2*σ_peak_2+p3*σ_peak_3))
    println(dKe*sum(σ_est_1+σ_est_2+σ_est_3),"\n")

    plot(be,p1*σ_peak_1+p2*σ_peak_2+p3*σ_peak_3)
    scatter(be,σ_est_1+σ_est_2+σ_est_3)

    plot(be,(dictAllData[symbol_h][1][:Snoisy]-dictAllData[symbol_h][1][:Sbg])/(dKe*sum(dictAllData[symbol_h][1][:Snoisy]-dictAllData[symbol_h][1][:Sbg])))
end

if false

    #TODO: remove noise estimation from this file
    ##
    ## svd noise estimation
    ##

    j_slect = 1
    λ_sym = Symbol(string("λe_",string(1.0e3λe[j_slect])));
    ν_sym = Symbol(string("hν_",Int64(round(hν[j_slect]))));
    F_λ  = svd(dictAllData[ν_sym][1][:σ_cs_dens]*dictAllGeom[λ_sym][1][:H]');

    bg_rm = 1.0
    s_th = 2

    UF_ν  = F_λ.U[:,s_th:end]'*(dictAllData[ν_sym][1][:Snoisy]-bg_rm*dictAllData[ν_sym][1][:Sbg]);

    noise_ν = F_λ.U[:,s_th:end]*UF_ν;

    σ_ν  = dictAllData[ν_sym][1][:Snoisy]-noise_ν;
    σ_ν[σ_ν.<=0.0] .= 0.0

    x_symbol = :Ke;
    figure(); 
    plot(dictAllData[ν_sym][1][x_symbol],dictAllData[ν_sym][1][:Stot]) 
    scatter(dictAllData[ν_sym][1][x_symbol],σ_ν)


    figure()
    scatter(collect(1:241),noise_ν,color="tab:blue")
    fill_between(collect(1:241),-sqrt.(σ_ν),sqrt.(σ_ν),alpha=0.5,color="tab:blue")


    # τ_ν = (σ_ν-dictAllData[ν_sym][1][:Sbg])./(dictAllData[ν_sym][1][:σ_cs_dens]*(dictAllGeom[λ_sym][1][:H]'*ρA_1))

    # figure()
    # plot(τ_ν[σ_ν.>0.1*maximum(σ_ν)])
    # ylim(0.0)

    # dictAllData[ν_sym][1][:σ_tot].*dictAllData[ν_sym][1][:F].*dictAllData[ν_sym][1][:T]/mean(τ_ν[σ_ν.>0.1*maximum(σ_ν)])

end


##
## alignement parameter estimation
##
# α_365,noise_ν = noiseAndParameterEstimation(dictAllData[:hν_365][1][:σ_cs_dens],dictAllGeom[Symbol("λe_0.5")][1][:H],Array{Cdouble,1}(dictAllData[:hν_365][1][:Snoisy]),dictAllData[:hν_365][1][:Sbg],dictAllGeom[Symbol("λe_0.5")][1][:ρ])
if false
    xc_vec = [-100.0; -75.0; -50.0; -20.0; -10.0; -5.0; -2.0; -1.0; 0.0; 1.0; 2.0; 5.0; 10.0; 20.0; 50.0; 75.0; 100.0; 125.0; 150.0; 175.0; 200.0; 250.0; 300.0; 400.0; 500.0];
    xc_vec = unique(sort([xc; μ0*sin(θ0); xc_vec]));

    α_al_off    = zeros(Cdouble,length(xc_vec));
    α_al_off_gt = zeros(Cdouble,length(xc_vec));
    for ic in 1:length(xc_vec)
        local bp = beamProfile(xc_vec[ic],yc,σx,σy)
        H_r,H_rθy,H_r_ph,H_rθy_ph,_,_,_,α_al_off[ic] =  alignmentParameter(bp,r,θ,y,x0,y0,z0,μ0,λe[1]);
        α_al_off_gt[ic] = (H_r_ph'*ρA_1)/(H_r'*ρA_1);
    end

    τ_al_noise    = zeros(Cdouble,Ndata);
    τ_al_noise_gt = zeros(Cdouble,Ndata);
    α_al_noise    = zeros(Cdouble,Ndata);
    α_al_noise_gt = zeros(Cdouble,Ndata);
    for i in 1:Ndata
        local symbol_h = Symbol(string("hν_",Int64(round(hν[i]))))
        local simbol_λ = Symbol(string("λe_",1.0e3λe[i]))
        τ_al_noise[i],_    = noiseAndParameterEstimation(dictAllData[symbol_h][1][:σ_cs_dens],dictAllGeom[simbol_λ][1][:H],Array{Cdouble,1}(dictAllData[symbol_h][1][:Snoisy]),dictAllData[symbol_h][1][:Sbg],ones(Cdouble,Nr))
        τ_al_noise_gt[i],_ = noiseAndParameterEstimation(dictAllData[symbol_h][1][:σ_cs_dens],dictAllGeom[simbol_λ][1][:H],Array{Cdouble,1}(dictAllData[symbol_h][1][:Snoisy]),dictAllData[symbol_h][1][:Sbg],ρA_1)
        α_al_noise[i] = τ_al_noise[i]/(Tj[i]*Fνj[i]*XPSpack.σ_C1s_interp[hν[i]])
        α_al_noise_gt[i] = τ_al_noise_gt[i]/(Tj[i]*Fνj[i]*XPSpack.σ_C1s_interp[hν[i]])
    end

    # τ_al_noise includes α_al_noise

    figure(); 
    plot(xc_vec,α_al_off_gt,color="tab:green")
    scatter(xc_vec,α_al_off,color="tab:blue")
    idx_sim = findfirst(xc_vec.==xc)
    idx_best = findfirst(xc_vec.==μ0*sin(θ0))
    scatter(μ0*sin(θ0),α_al_off_gt[idx_best],color="tab:red")
    scatter(xc*ones(Cdouble,Ndata),α_al_noise_gt,color="tab:olive")
    scatter(xc*ones(Cdouble,Ndata),α_al_noise,color="tab:orange")
    xlabel("horizontal deviation [\$\\mu m\$]",fontsize=14)
    ylabel("alignment [\$\\mu m^{-2}\$]",fontsize=14)
    xticks(fontsize=12)
    yticks(fontsize=12)
    ticklabel_format(axis="y",style="sci",scilimits=(-2,2))
    legend(["approximation","GT","closest to analyzer","from data with GT","from data"])
    tight_layout(pad=1.0, w_pad=0.5, h_pad=0.2)

    # savefig("alignment_factor.png")
    # savefig("alignment_factor.pdf")

end