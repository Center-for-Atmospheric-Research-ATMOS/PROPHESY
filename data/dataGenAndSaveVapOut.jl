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

FLAG_0001 = true               # selection of the profile (one must be true and the others false)
FLAG_0002 = false
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
hν = collect(LinRange(hν1, hν2,Ndata));                # photon energy range

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
# j = Ndata; # 1; # 10;                                               # select the photon energy
θ_aperture = 0.5*π/4
α_Ω = 4π*sin(θ_aperture/2.0)^2
# hν = collect(LinRange(hν1, hν2,Ndata));                         # central photon energy for each measurement
dhν = hν.*((1.0/25000.0)*(hν.<500.0) + (1.0/15000.0)*(hν.>=500.0)); # bandwidth of the photon beam
Fνj = 1.5e11*ones(Cdouble,Ndata); #  1.0e3*ones(Cdouble,Ndata); # 3.0e11*α_al; #                                     # flux densities
Tj   = α_Ω*(10.0.+0.0collect(LinRange(5.0,10.0,Ndata))); # LinRange(5.0,10.0,Ndata))                          # transmission factors
σ_ke = 2.0*dKe*ones(Cdouble,Ndata); # collect(LinRange(1.0,2.0,Ndata)); #                      # kinetic energy bandwidths of the analyzer (one per photon energy)


# dictionary where to push the data and geometry factor
dictAllData = Dict()
dictAllGeom = Dict()
for j in 1:Ndata # can potentially multi-thread this loop: but need to sync before writing files
    # for a given photon energy νj, measure a spectrum
    local Keij = reverse(hν[j] .- Be) ;                                       # centers of the analyzer's channels
    local kKe = floor(5σ_ke[j]/dKe);                                             # number of extra discretization point so that the channels at Ki[1] and Ki[end] are not missing to much information
    local Kdisc = [Keij[1].-dKe*reverse(collect(1:kKe)); Keij; Keij[end].+dKe*collect(1:kKe)]; # discretization point
    local μKe = 0.5*(Keij[1]+Keij[end]);                                      # central kinetic energy (a bit the same role as pass energy)
    # local σ_bg = 0.05*(hν[j]-BeC1s) # μKe;
    ##
    ## compute the cross section of the sample (to be estimated from the data in an estimation setting)
    ##

    local Be0 = hν[j] .- Keij;
    local σ_cs_fg =  σ_cs.(hν[j],Keij,μKe); # ./XPSpack.σ_C1s_interp[hν[j]] 

    ##
    ## electron flux without the geometry factor: signal of interest (Sj) and background signal (Sj_bg)
    ##
    local Sj,Sj_bg,_,_ = simulateSpectrum(Fνj[j],hν[j],dhν[j],
        Keij,σ_ke[j],Tj[j],Kdisc,
        reverse(Be0),reverse(σ_cs_fg),σ_bg(μKe)*σ_bg_density(Keij,hν[j]-BeC1s,ΔBeC1s));

    ##
    ## geometry factors
    ##
    H_deom,_,H_geom,_,_,_,_,_ = alignmentParameter(bp,r,θ,y,x0,y0,z0,μ0,λe[j])

    ##
    ## for the signal without noise
    ##
    SbgC1s      = (H_geom'*ρA_1)*Sj_bg
    SpectrumA_1 = (H_geom'*ρA_1)*Sj

    ##
    ## add noise
    ##
    local SC1snoise = countElectrons(SbgC1s+SpectrumA_1)

    # plot signals w.r.t. the kinetic energy
    figure(); plot(Keij,SbgC1s); plot(Keij,SbgC1s+SpectrumA_1); scatter(Keij,rand.(Poisson.(SbgC1s+SpectrumA_1))); xlabel("kinetic energy [eV]"); ylabel("spectrum [a.u.]") 

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


if false

    #TODO: remove noise estimation from this file
    ##
    ## svd noise estimation
    ##


    F_365  = svd(dictAllData[:hν_365][1][:σ_cs_dens]*dictAllGeom[Symbol("λe_0.5")][1][:H]')
    F_649  = svd(dictAllData[:hν_649][1][:σ_cs_dens]*dictAllGeom[Symbol("λe_1.75")][1][:H]')
    F_932  = svd(dictAllData[:hν_932][1][:σ_cs_dens]*dictAllGeom[Symbol("λe_3.0")][1][:H]')
    F_1216 = svd(dictAllData[:hν_1216][1][:σ_cs_dens]*dictAllGeom[Symbol("λe_4.25")][1][:H]')
    F_1500 = svd(dictAllData[:hν_1500][1][:σ_cs_dens]*dictAllGeom[Symbol("λe_5.5")][1][:H]')


    F_365.U'*(dictAllData[:hν_365][1][:Snoisy]-dictAllData[:hν_365][1][:Sbg])


    bg_rm = 1.0
    s_th = 2

    UF_365  = F_365.U[:,s_th:end]'*(dictAllData[:hν_365][1][:Snoisy]-bg_rm*dictAllData[:hν_365][1][:Sbg]);
    UF_649  = F_649.U[:,s_th:end]'*(dictAllData[:hν_649][1][:Snoisy]-bg_rm*dictAllData[:hν_649][1][:Sbg]);
    UF_932  = F_932.U[:,s_th:end]'*(dictAllData[:hν_932][1][:Snoisy]-bg_rm*dictAllData[:hν_932][1][:Sbg]);
    UF_1216 = F_1216.U[:,s_th:end]'*(dictAllData[:hν_1216][1][:Snoisy]-bg_rm*dictAllData[:hν_1216][1][:Sbg]);
    UF_1500 = F_1500.U[:,s_th:end]'*(dictAllData[:hν_1500][1][:Snoisy]-bg_rm*dictAllData[:hν_1500][1][:Sbg]);

    noise_365  = F_365.U[:,s_th:end]*UF_365
    noise_649  = F_649.U[:,s_th:end]*UF_649
    noise_932  = F_932.U[:,s_th:end]*UF_932
    noise_1216 = F_1216.U[:,s_th:end]*UF_1216
    noise_1500 = F_1500.U[:,s_th:end]*UF_1500

    σ_365  = dictAllData[:hν_365][1][:Snoisy]-noise_365;
    σ_365[σ_365.<=0.0] .= 0.0
    σ_649  = dictAllData[:hν_649][1][:Snoisy]-noise_649; 
    σ_649[σ_649.<=0.0] .= 0.0
    σ_932  = dictAllData[:hν_932][1][:Snoisy]-noise_932;
    σ_932[σ_932.<=0.0] .= 0.0
    σ_1216 = dictAllData[:hν_1216][1][:Snoisy]-noise_1216; 
    σ_1216[σ_1216.<=0.0] .= 0.0
    σ_1500 = dictAllData[:hν_1500][1][:Snoisy]-noise_1500; 
    σ_1500[σ_1500.<=0.0] .= 0.0

    x_symbol = :Ke;
    figure(); 
    plot(dictAllData[:hν_365][1][x_symbol],dictAllData[:hν_365][1][:Stot]) # collect(1:241)
    scatter(dictAllData[:hν_365][1][x_symbol],σ_365)

    plot(dictAllData[:hν_649][1][x_symbol],dictAllData[:hν_649][1][:Stot])
    scatter(dictAllData[:hν_649][1][x_symbol],σ_649)

    plot(dictAllData[:hν_932][1][x_symbol],dictAllData[:hν_932][1][:Stot])
    scatter(dictAllData[:hν_932][1][x_symbol],σ_932)

    plot(dictAllData[:hν_1216][1][x_symbol],dictAllData[:hν_1216][1][:Stot])
    scatter(dictAllData[:hν_1216][1][x_symbol],σ_1216)

    plot(dictAllData[:hν_1500][1][x_symbol],dictAllData[:hν_1500][1][:Stot])
    scatter(dictAllData[:hν_1500][1][x_symbol],σ_1500)



    figure()
    scatter(collect(1:241),noise_365,color="tab:blue")
    fill_between(collect(1:241),-sqrt.(σ_365),sqrt.(σ_365),alpha=0.5,color="tab:blue")

    scatter(collect(1:241),noise_649,color="tab:orange")
    fill_between(collect(1:241),-sqrt.(σ_649),sqrt.(σ_649),alpha=0.5,color="tab:orange")

    scatter(collect(1:241),noise_932,color="tab:green")
    fill_between(collect(1:241),-sqrt.(σ_932),sqrt.(σ_932),alpha=0.5,color="tab:green")

    scatter(collect(1:241),noise_1216,color="tab:red")
    fill_between(collect(1:241),-sqrt.(σ_1216),sqrt.(σ_1216),alpha=0.5,color="tab:red")

    scatter(collect(1:241),noise_1500,color="tab:purple")
    fill_between(collect(1:241),-sqrt.(σ_1500),sqrt.(σ_1500),alpha=0.5,color="tab:purple")

    # τ_365 = (σ_365-dictAllData[:hν_365][1][:Sbg])./(dictAllData[:hν_365][1][:σ_cs_dens]*sum(dictAllGeom[Symbol("λe_0.5")][1][:H][1:183]))
    τ_365 = (σ_365-dictAllData[:hν_365][1][:Sbg])./(dictAllData[:hν_365][1][:σ_cs_dens]*(dictAllGeom[Symbol("λe_0.5")][1][:H]'*ρA_1))

    figure()
    plot(τ_365[σ_365.>0.1*maximum(σ_365)])
    ylim(0.0)

    dictAllData[:hν_365][1][:σ_tot].*dictAllData[:hν_365][1][:F].*dictAllData[:hν_365][1][:T]/mean(τ_365[σ_365.>0.1*maximum(σ_365)])

end


##
## alignement parameter estimation
##
# α_365,noise_365 = noiseAndParameterEstimation(dictAllData[:hν_365][1][:σ_cs_dens],dictAllGeom[Symbol("λe_0.5")][1][:H],Array{Cdouble,1}(dictAllData[:hν_365][1][:Snoisy]),dictAllData[:hν_365][1][:Sbg],dictAllGeom[Symbol("λe_0.5")][1][:ρ])
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

# # TODO: check again methodology here... it does not seem to be working properly or at least not as intended
# function parameterEstimation(σ_χ::Array{Cdouble,1},H::Array{Cdouble,1},I_data::Array{Cdouble,1},I_bg::Array{Cdouble,1},ρ::Array{Cdouble,1})
#     F  = svd(σ_χ*H')
#     noise_data = F.U[:,2:end]*(F.U[:,2:end]'*(I_data-I_bg));
#     σ_data = I_data-(I_bg+noise_data)

#     # figure(); plot(I_data); plot(I_data-noise_data); plot(σ_data)
 
#     mean(σ_data[σ_χ.>0.1*maximum(σ_χ)]./(σ_χ[σ_χ.>0.1*maximum(σ_χ)]*(H'*ρ))),noise_data
#  end

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