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
SHORT_RANGE = true # WARNING: not ready for wide range (i.e. SHORT_RANGE=false)

MODEL_5   = false               # select the number of attenuation lengths probed
MODEL_10  = false
MODEL_20  = true

FLAG_OFF_CENTER_0 = true;
FLAG_OFF_CENTER_1 = false;
FLAG_OFF_CENTER_2 = false;
FLAG_OFF_CENTER_3 = false;
FLAG_OFF_CENTER_4 = false;
FLAG_OFF_CENTER_5 = false;

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
x0 = sqrt(2.0)*5000.0 # (x0,y0,z0) are the coordinates of the analyzer's apperture
y0 = 0.0;
z0 = 5000.0;

# soon add the option of sphere or plane geometry?
save_folder = string(save_folder,"cylinder_radius_",μ0,"/peak_shift/")


# spacial discretization 
r = collect(range(μ0-k0*λe0,μ0+δr,length=Nr));
r_lowres = collect(range(μ0-k0*λe0,μ0+δr,length=Nr_lowres));
θ0 = atan(x0,z0)
θ = collect(range(θ0-π/2.0,θ0+π/2.0,Nθ));
y = collect(range(-L/2.0,L/2.0,length=Ny));

r_gas = collect(range(μ0,10μ0,length=Nr));
θ_gas = collect(range(θ0-π/2.0-acos(μ0/r_gas[end]),θ0+π/2.0+acos(μ0/r_gas[end]),Nθ));
y_gas = y;

# concentration profiles (4 different cases)
ρ0 = 1.0
ρ_vac = 0.0
ρ_bulk = 55.49; # M (water concentration in liquid phase)

M_water = 18.02e-3 # kg/mol
R_gas = 8.314   # J/(K*mol)
T_gas = 280.0   # K
P_gas = (0.5/760.0)*1.0e5   # Pa (0.5 Torr?)
ρ_gas = (P_gas*M_water/(R_gas*T_gas))*(μ0./r_gas); # assume that the gas concentration profile is decreasing inversly proportional to the distance to the sample (the amount of mater crossing the surface of the cylinder is the same that crossing a bigger virtual cylinder)


ρA_1 = ρ_bulk*logistic.(1000.0*(μ0.-r)/0.5,ρ_vac,ρ0,1.0);
exp_tag     = "water"

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

# not sure it's that short anymore
if SHORT_RANGE                                          # bounds of the attenuation lengths
    λe1 = 1.3; hν1 =  650.0;
    λe2 = 3.8; hν2 = 1884.0;
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
σx = 100.0; # spread of the beam
σy = 25.0;
xc = 100.0; # center of the beam
yc = 99.0;
if FLAG_OFF_CENTER_0
    kc=0
end
if FLAG_OFF_CENTER_1
    kc=1
end
if FLAG_OFF_CENTER_2
    kc=2
end
if FLAG_OFF_CENTER_3
    kc=3
end
if FLAG_OFF_CENTER_4
    kc=4
end
if FLAG_OFF_CENTER_5
    kc=5
end
xc = kc*σx;
yc = 99.0;
save_folder = string(save_folder,"offcenter_",kc,"/")
# beam profile
bp = beamProfile(xc,yc,σx,σy);


## 
## cross section density model: just a dummy model looking like O1s
## 
dKe = 0.05;
Be = collect(528.0:dKe:542.0);
Nspectrum = length(Be);
BeO1s = mean(Be);
ΔBeO1s = 3.0;

hknot_peak = [650.0;  651.0;  789.0;      897.0; 907.0;      1154.0;   1197.0; 1565.0];
Eb_knot_liq = [533.84;  533.882; 534.392; 534.456; 534.581; 534.869; 534.667; 533.989];
Eb_knot_gas = [535.818; 535.756; 536.273; 536.324; 536.438; 536.680; 536.479; 535.765];

σEb_knot_liq = 0.5*[1.091; 1.008 ; 1.498; 1.608; 1.618; 1.624; 1.659; 1.641];
σEb_knot_gas = 0.5*[0.629; 0.631; 0.804; 0.853;  0.863; 0.949; 0.995; 1.016];

Ahν = [hknot_peak.^3 hknot_peak.^2 hknot_peak.^1 hknot_peak.^0];

μBe_var_liq = (inv(Ahν'*Ahν)*Ahν'*Eb_knot_liq)';
μBe_var_gas = (inv(Ahν'*Ahν)*Ahν'*Eb_knot_gas)';

σBe_var_liq = (inv(Ahν'*Ahν)*Ahν'*σEb_knot_liq)';
σBe_var_gas = (inv(Ahν'*Ahν)*Ahν'*σEb_knot_gas)';


hhν = collect(310.0:20.0:1900.0);
figure()
scatter(hknot_peak,Eb_knot_liq)
plot(hhν,dropdims(μBe_var_liq*[hhν.^3 hhν.^2 hhν.^1 hhν.^0]',dims=1))
scatter(hknot_peak,Eb_knot_gas)
plot(hhν,dropdims(μBe_var_gas*[hhν.^3 hhν.^2 hhν.^1 hhν.^0]',dims=1))

figure()
scatter(hknot_peak,σEb_knot_liq)
plot(hhν,dropdims(σBe_var_liq*[hhν.^3 hhν.^2 hhν.^1 hhν.^0]',dims=1))
scatter(hknot_peak,σEb_knot_gas)
plot(hhν,dropdims(σBe_var_gas*[hhν.^3 hhν.^2 hhν.^1 hhν.^0]',dims=1))

function σ_cs_O1s_liq(hν::Cdouble,Ke::Array{Cdouble,1})
    Be = hν.-Ke;
    # make the peaks vary in location and width
    μBe_O1s_liq = μBe_var_liq*[hν^3; hν^2; hν; 1.0];
    σBe_O1s_liq = σBe_var_liq*[hν^3; hν^2; hν; 1.0];
    # partial cross section (one for each chemical state)
    σ_peak = (1.0/sqrt(2.0π*σBe_O1s_liq^2))*exp.(-(Be.-μBe_O1s_liq).^2/(2.0σBe_O1s_liq^2));

    # cross section value
    XPSpack.σ_O1s_interp[hν]*σ_peak,μBe_O1s_liq,σBe_O1s_liq
end

function σ_cs_O1s_gas(hν::Cdouble,Ke::Array{Cdouble,1})
    Be = hν.-Ke;
    # make the peaks vary in location and width
    μBe_O1s_gas = μBe_var_gas*[hν^3; hν^2; hν; 1.0];
    σBe_O1s_gas = σBe_var_gas*[hν^3; hν^2; hν; 1.0];
    # partial cross section (one for each chemical state)
    σ_peak = (1.0/sqrt(2.0π*σBe_O1s_gas^2))*exp.(-(Be.-μBe_O1s_gas).^2/(2.0σBe_O1s_gas^2));

    # cross section value
    XPSpack.σ_O1s_interp[hν]*σ_peak,μBe_O1s_gas,σBe_O1s_gas
end

##
## acqusition parameters
##

θ_aperture = 0.5*π/4                                                       # aperture's cone angle
α_Ω        = 4π*sin(θ_aperture/2.0)^2;                                     # aperture's solid angle
Tj         = α_Ω*ones(Cdouble,Ndata);                                      # transmission factors
σ_ke       = 2.0*dKe*ones(Cdouble,Ndata);                                  # kinetic energy bandwidths of the analyzer (one per photon energy)
dhν        = hν.*((1.0/25000.0)*(hν.<500.0) + (1.0/15000.0)*(hν.>=500.0)); # bandwidth of the photon beam

hknot = [650.0; 1884.0; 1315.0; 907.0]
Fknot = [3.943312e+13; 1.419204e+14; 3.853618e+14; 1.651883e+14]
μF_var =  (inv([hknot.^3 hknot.^2 hknot hknot.^0])*Fknot)';
Fνj    = dropdims(μF_var*[hν.^3 hν.^2 hν.^1 hν.^0]',dims=1);
# figure()
# scatter(hν,Fνj)
# scatter(hknot,Fknot)
# scatter(collect(300.0:10.0:2000.0),dropdims(μF_var*[collect(300.0:10.0:2000.0).^3 collect(300.0:10.0:2000.0).^2 collect(300.0:10.0:2000.0).^1 collect(300.0:10.0:2000.0).^0]',dims=1))



# dictionary where to push the data and geometry factor
dictAllData  = Dict()
dictAllGeom  = Dict()
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
    local σ_cs_liq,μBe_liq,σ_be_liq = σ_cs_O1s_liq(hν[j],Keij)
    local σ_cs_gas,μBe_gas,σ_be_gas = σ_cs_O1s_gas(hν[j],Keij)

    ##
    ## electron flux without the geometry factor: signal of interest (Sj) and background signal (Sj_bg)
    ##
    local Sj,Sj_bg,_,_ = simulateSpectrum(Fνj[j],hν[j],dhν[j],
        Keij,σ_ke[j],Tj[j],Kdisc,
        reverse(Be0),reverse(σ_cs_liq),σ_bg(μKe)*σ_bg_density(Keij,hν[j]-BeO1s,ΔBeO1s));

    ##
    ## geometry factors
    ##
    local H_deom,_,H_geom,_,_,_,_,al_liq = alignmentParameter(bp,r,θ,y,x0,y0,z0,μ0,λe[j])

    ##
    ## for the signal without noise
    ##
    local SbgO1s      = (H_geom'*ρA_1)*Sj_bg
    local SpectrumA_1 = (H_geom'*ρA_1)*Sj

    ##
    ## gas phase signal
    ##
    local H_deom_gas,_,H_geom_gas,_,_,_,_,al_gas = alignmentParameter(bp,r_gas,θ_gas,y_gas,x0,y0,z0,μ0,λe[j])
    
    local Sj_gas,_,_,_ = simulateSpectrum(Fνj[j],hν[j],dhν[j],
        Keij,σ_ke[j],Tj[j],Kdisc,
        reverse(Be0),reverse(σ_cs_gas),zeros(Cdouble,Nspectrum));
    local SpectrumA_1_gas = (H_geom_gas'*ρ_gas)*Sj_gas;

    # r_gas = collect(range(μ0,5μ0,length=Nr));
    # θ_gas = collect(range(θ0-π/2.0-acos(μ0/r_gas[end]),θ0+π/2.0+acos(μ0/r_gas[end]),Nθ));
    # y_gas = y;
    # H_r,H_rθy,H_r_ph,H_rθy_ph,_,_,_,_ = alignmentParameter(bp,r_gas,θ_gas,y_gas,x0,y0,z0,μ0,λe[1])
    # figure(); imshow(H_rθy_ph[:,:,128]); colorbar()
    # ax1 = subplot(121,polar=true)
    # pcm1 = ax1.pcolormesh(θ_gas,r_gas,H_rθy_ph[:,:,128].*(0.5μ0./r_gas),edgecolors="face")
    # ax2 = subplot(122,polar=true)
    # pcm2 = ax2.pcolormesh(θ_gas,r_gas,H_rθy_ph[:,:,1].*(0.5μ0./r_gas),edgecolors="face")

    ##
    ## add noise
    ##
    local SO1snoise = countElectrons(SbgO1s+SpectrumA_1+SpectrumA_1_gas)


    ##
    ## push data to dicts
    ##
    
    local dictData = Dict( "Ke" => Keij, "Be" => Be0, "μKe" => μKe,
        "σ_cs_dens" => σ_cs_liq./XPSpack.σ_O1s_interp[hν[j]], "σ_cs_dens_gas" => σ_cs_gas./XPSpack.σ_O1s_interp[hν[j]], "σ_tot" => XPSpack.σ_O1s_interp[hν[j]],
        "SpectrumA_1" => SpectrumA_1, "SpectrumA_1_gas" => SpectrumA_1_gas, "Sbg" => SbgO1s, 
        "Stot" => SbgO1s+SpectrumA_1, "Snoisy" => SO1snoise,
        "T" => Tj[j], "λ" => 1.0e3λe[j], "F" => Fνj[j], "hν" => hν[j], "α" => sum(SpectrumA_1)/sum(SpectrumA_1_gas));
    local dictMetaData = Dict("μKe" => μKe,
        "σ_tot" => XPSpack.σ_O1s_interp[hν[j]], "peak_mode_liq" => μBe_liq, "peak_width_liq"=>σ_be_liq, "peak_mode_gas" => μBe_gas, "peak_width_gas"=>σ_be_gas,
        "T" => Tj[j], "λ" => 1.0e3λe[j], "F" => Fνj[j], "hν" => hν[j]);

    local dictGeom = Dict("model" => "sharp edge cylinder + outside vapor", 
                "hν" => hν[j], "λ" => 1.0e3λe[j], "radius" => μ0, "max_depth" => k0*λe0,
                    "x0" => x0, "y0" => y0, "z0" => z0, "δr" => δr, "xc" => xc, "yc" => yc, "σx" => σx, "σy" => σy, "α" => al_liq,
                    "r" => r, "r_gas" => r_gas, "H" => H_deom, "H_true" => H_geom, "H_gas"=>H_deom_gas, "H_gas_true"=>H_geom_gas, 
                    "ρ" => ρA_1, "ρ_gas" => ρ_gas)

    local dfGeom     = DataFrame(dictGeom);
    local dfData     = DataFrame(dictData);
    local dfMetaData = DataFrame(dictMetaData);

    dictAllData[Symbol(string("hν_",Int64(round(hν[j]))))]      = (eachcol(dfData),names(dfData))
    dictAllData[Symbol(string("meta_hν_",Int64(round(hν[j]))))] = (eachcol(dfMetaData),names(dfMetaData))
    dictAllGeom[Symbol(string("λe_",string(1.0e3λe[j])))]       = (eachcol(dfGeom),names(dfGeom))
end

if SAVE_DATA
    mkpath(string(save_folder,exp_tag,"/"));
    XLSX.writetable(string(save_folder,exp_tag,"/data.xlsx"); dictAllData...) 
    XLSX.writetable(string(save_folder,exp_tag,"/model.xlsx"); dictAllGeom...)
end
include("plotData.jl")

