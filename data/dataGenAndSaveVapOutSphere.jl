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

FLAG_OFF_CENTER_0 = false;
FLAG_OFF_CENTER_1 = false;
FLAG_OFF_CENTER_2 = false;
FLAG_OFF_CENTER_3 = true;

FLAG_0001 = false               # selection of the profile (one must be true and the others false)
FLAG_0002 = false
FLAG_0003 = false
FLAG_0004 = true

save_folder = "./";
SAVE_DATA = false   # flag for saving data and model
SAVE_FIG  = false   # save the simulated data or not
SWEEPS_ON = false
nb_sweeps = 50;

# geometry setup
λe0 = 2.0e-3;        # reference penetration depth in μm
δr = 10.0e-3          # transition to vacuum layer thickness (let's set about 1 nm)
k0 = 10;             # compute the measurement model over a distance of k0*λe0
Nr = 201;            # number of discretization points in the radial dimension
Nφ = 256;            # number of discretization points in the polar angle dimension
Nθ = 257;            # number of discretization points in the azimuthal angle dimension
μ0 = 10.0; # 20.0;   # radius of the sphere

# L = 50.0 # 60.0 # 20; # 100.0;     # height of the irradiated sample (the vertical extent of the beam is more like 20μm instead of 100μm)
x0 = sqrt(2.0)*5000.0 # (x0,y0,z0) are the coordinates of the analyzer's apperture
y0 = 0.0;
z0 = 5000.0;

# soon add the option of sphere or plane geometry?
save_folder = string(save_folder,"sphere_radius_",μ0,"/peak_shift/")


# spacial discretization 
r = collect(range(μ0-k0*λe0,μ0+δr,length=Nr));
# θ0 = atan(x0,z0);
# φ0 = acos(y0/sqrt(x0^2 + y0^2 + z0^2));
# φ = collect(range(φ0-π/2.0,φ0+π/2.0,Nφ)); 
# θ = collect(range(θ0-π/2.0,θ0+π/2.0,Nθ));
φ = collect(range(0.0,π,Nφ)); 
θ = collect(range(0.0,2.0π,Nθ));

# concentration profiles (4 different cases)
ρ0 = 1.0
ρ_vac = 0.0

# ρ_bulk =  5.0e-3; #  5mM
ρ_bulk = 10.0e-3; # 10mM


if FLAG_0001
    ρA_1 = ρ_bulk*logistic.(1000.0*(μ0.-r)/2.5,ρ_vac,ρ0,1.0);
    exp_tag     = "0001"
end
if FLAG_0002
    ρA_1 = ρ_bulk*(logistic.(1000.0*(μ0.-r)/2.5,ρ_vac,ρ0,1.0) .+ 2.0exp.(-(1000.0*(μ0.-r).-0.0).^2. /(2.0*1.25^2)));
    exp_tag     = "0002"
end
if FLAG_0003
    ρA_1 = ρ_bulk*(logistic.(1000.0*(μ0.-r)/2.5,ρ_vac,ρ0,1.0) .+ exp.(-(1000.0*(μ0.-r).-0.5).^2. /(2.0*2.5^2)));
    exp_tag     = "0003"
end
if FLAG_0004
    ρA_1 = ρ_bulk*exp.(-(1000.0((μ0.-r))).^2. /(2.0*2.5^2));
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

# not sure it's that short anymore
if SHORT_RANGE                                          # bounds of the attenuation lengths
    # λe1 = 1.3; hν1 = 365.0;
    # λe2 = 2.5; hν2 = 1500.0;
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
σx = 0.5*100.0; # spread of the beam
σy = 0.5*25.0;
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
xc = kc*σx;
yc = 99.0;

xc = kc*σx;
yc = 50.0 # 40.0 # 75.0; # 5.0*σy; # 6.0*σy #WARNING: this value accounts for the 

save_folder = string(save_folder,"offcenter_",kc,"/")
# beam profile
bp = beamProfile(xc,yc,σx,σy);
α_al = zeros(Cdouble,Ndata);
for i in 1:Ndata
    local H_r,_,H_r_ph,_,_,_,_,α_al[i] =  alignmentParameterSphere(bp,r,φ,θ,x0,y0,z0,μ0,λe[i]);
end

## 
## cross section density model: just a dummy model looking like C1s
## 

dKe = 0.05;
# Be = collect(286.0:dKe:298.0);
Be = collect(271.0:dKe:291.0); #
Nspectrum = length(Be);
# μBe = [290.2; 292.0; 293.0]
# μBe = [280.436245; 281.051876; 283.813405]
μBe = [280.3; 281.8; 283.8]

# σ_be = sqrt(2.0)*[0.45; 0.25; 0.6];
σ_be = 0.5e-3*[935.238187; 935.238187; 854.723209];

BeC1s = 284.0 # mean(Be);
ΔBeC1s = 6.0 # 12.0/5.0#(Be[end]-Be[1])/5;

# The width may change, let's say, from 1000 to 1300 a.u on depth range 1,5...3,7 nm, for example. Peak positions may vary up to 4-5 nm, in my case it's from 284 to 280 eV for a SDS head group on depth range mentioned above. 

Ahν  = [650.0^2  650.0 1.0; 
1315.0^2 1315.0 1.0;
1884.0^2 1884.0 1.0];
Bhν = [1.0; (286.0+1.0)/286.0; (286.0-3.0)/286.0];
μBe_var = inv(Ahν)*Bhν;

hhν = collect(310.0:20.0:1900.0);
figure()
scatter([650.0; 1315.0; 1884.0],Bhν)
plot(hhν,dropdims(μBe_var'*[hhν.^2 hhν.^1 hhν.^0]',dims=1))


# figure();
# plot(collect(650.0:10.0:1884.0),dropdims(μBe_var'*[collect(650.0:10.0:1884.0).^2 collect(650.0:10.0:1884.0) ones(Cdouble,length(collect(650.0:10.0:1884.0)))]',dims=1))

function σ_cs(hν::Cdouble,Ke::Array{Cdouble,1},μKe::Cdouble;μKe0::Cdouble=50.0,μKe1::Cdouble=1200.0)
    Be = hν.-Ke;
    # make the peaks vary in location and width
    peak_dev = μBe_var'*[hν^2; hν; 1.0];
    μBe_dev = peak_dev*μBe
    # σ_be_var = (1.0+0.3*((μKe-μKe0)/(μKe1-μKe0)))*σ_be;
    σ_be_var = (1.0+0.3*((hν-650.0)/(1884.0-650.0)))*σ_be;
    # partial cross section (one for each chemical state)
    σ_peak_1 = (1.0/sqrt(2.0π*σ_be_var[1]^2))*exp.(-(Be.-μBe_dev[1]).^2/(2.0σ_be_var[1]^2));
    σ_peak_2 = (1.0/sqrt(2.0π*σ_be_var[2]^2))*exp.(-(Be.-μBe_dev[2]).^2/(2.0σ_be_var[2]^2));
    σ_peak_3 = (1.0/sqrt(2.0π*σ_be_var[3]^2))*exp.(-(Be.-μBe_dev[3]).^2/(2.0σ_be_var[3]^2));
    # quantity of chemical states
    p1 = 0.7 # 0.85 .+ (0.77-0.85)*(μKe.-μKe0)./(μKe1-μKe0); # 0.8 # 
    p2 = 0.25 # 0.125 .+ (0.12-0.125)*(μKe.-μKe0)./(μKe1-μKe0); # 0.12 # 
    p3 = 1.0-(p1+p2);

    # cross section value (only for hν ∈ [295,1500.0])
    XPSpack.σ_C1s_interp[hν]*(p1*σ_peak_1+p2*σ_peak_2+p3*σ_peak_3),μBe_dev,σ_be_var,[p1;p2;p3]
end

##
## acqusition parameters
##

θ_aperture = 0.5*π/4                                                       # aperture's cone angle
α_Ω        = 4π*sin(θ_aperture/2.0)^2;                                     # aperture's solid angle
# Tj         = α_Ω*(10.0.+0.0collect(LinRange(5.0,10.0,Ndata)));             # transmission factors
Tj         = α_Ω*ones(Cdouble,Ndata);                                      # transmission factors
σ_ke       = 2.0*dKe*ones(Cdouble,Ndata);                                  # kinetic energy bandwidths of the analyzer (one per photon energy)
dhν        = hν.*((1.0/25000.0)*(hν.<500.0) + (1.0/15000.0)*(hν.>=500.0)); # bandwidth of the photon beam
# Fνj        = 1.5e11*ones(Cdouble,Ndata);                                   # flux densities
# Fνj        = 1.5e12*ones(Cdouble,Ndata);                                   # flux densities

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
    # local σ_cs_fg,μBe_dev,σ_be_var,τ_be =  σ_cs.(hν[j],Keij,μKe);
    local σ_cs_fg,μBe_dev,σ_be_var,τ_be =  σ_cs(hν[j],Keij,μKe);

    ##
    ## electron flux without the geometry factor: signal of interest (Sj) and background signal (Sj_bg)
    ##
    local Sj,Sj_bg,_,_ = simulateSpectrum(Fνj[j],hν[j],dhν[j],
        Keij,σ_ke[j],Tj[j],Kdisc,
        reverse(Be0),reverse(σ_cs_fg),σ_bg(μKe)*σ_bg_density(Keij,hν[j]-BeC1s,ΔBeC1s)); # σ_bg_lin_density(Keij,hν[j]-BeC1s,5.0e4ΔBeC1s) 5.0

    ##
    ## geometry factors
    ##
    local H_deom,_,H_geom,_,_,_,_,al_C1s = alignmentParameterSphere(bp,r,φ,θ,x0,y0,z0,μ0,λe[j])

    ##
    ## for the signal without noise
    ##
    local SbgC1s      = (H_geom'*ρA_1)*Sj_bg
    local SpectrumA_1 = (H_geom'*ρA_1)*Sj

    ##
    ## add noise
    ##
    local SC1snoise = countElectrons(SbgC1s+SpectrumA_1,1.0)

    if SWEEPS_ON # just repeat the acquisition
        local SC1snoise_sweeps = zeros(Cdouble,length(SC1snoise),nb_sweeps)
        for i in 1:nb_sweeps
            SC1snoise_sweeps[:,i] = countElectrons(SbgC1s+SpectrumA_1,0.0);
        end
    end

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
    if SWEEPS_ON
        [dictData[string("Snoisy_",i)] = SC1snoise_sweeps[:,i] for i in 1:nb_sweeps]
    end
    local dictMetaData = Dict("μKe" => μKe,
        "σ_tot" => XPSpack.σ_C1s_interp[hν[j]], "peak_mode" => μBe_dev, "peak_width"=>σ_be_var, "peak_probability"=>τ_be,
        "T" => Tj[j], "λ" => 1.0e3λe[j], "F" => Fνj[j], "hν" => hν[j]);

    local dictGeom = Dict("model" => "sharp edge cylinder + outside vapor", 
                "hν" => hν[j], "λ" => 1.0e3λe[j], "radius" => μ0, "max_depth" => k0*λe0,
                    "x0" => x0, "y0" => y0, "z0" => z0, "δr" => δr,"xc" => xc, "yc" => yc, "σx" => σx, "σy" => σy, "α" => al_C1s,
                    "r" => r, "H" => H_deom, "H_true" => H_geom, "ρ" => ρA_1)

    local dfGeom     = DataFrame(dictGeom);
    local dfData     = DataFrame(dictData);
    local dfMetaData = DataFrame(dictMetaData);

    dictAllData[Symbol(string("hν_",Int64(round(hν[j]))))]      = (eachcol(dfData),names(dfData))
    dictAllData[Symbol(string("meta_hν_",Int64(round(hν[j]))))] = (eachcol(dfMetaData),names(dfMetaData))
    dictAllGeom[Symbol(string("λe_",string(1.0e3λe[j])))]       = (eachcol(dfGeom),names(dfGeom))
end

if SAVE_DATA
    mkpath(string(save_folder,exp_tag,"/"));
    XLSX.writetable(string(save_folder,exp_tag,"/data.xlsx"); dictAllData...) # TODO: get outside the loop
    XLSX.writetable(string(save_folder,exp_tag,"/model.xlsx"); dictAllGeom...)
end
plot_sym=:Be
include("plotData.jl")
