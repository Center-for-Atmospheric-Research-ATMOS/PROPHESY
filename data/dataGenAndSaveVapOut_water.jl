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

# modeling XPS
using XPSpack # experiment model (geometry factor and cross section estimation)
using ATTIRE  # kinetic energy analyzer



# flags
MODEL_5   = true               # select the number of attenuation lengths probed
MODEL_10  = false
MODEL_20  = false


SAVE_DATA = false  # flag for saving data and model
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
L = 20.0 # 25; # 100.0;     # height of the irradiated sample (the vertical extent of the beam is more like 20μm instead of 100μm)
x0 = sqrt(2.0)*5000.0 # (x0,y0,z0) are the coordinates of the analyzer's apperture
y0 = 0.0;
z0 = 5000.0;

# integration inteval
Δt = 10.0 # [s]
D_dark = 0.0 # dark current coefficient is a number of electron per unit time [# s^{-1}]

# unit conversion constant (some of the quantities are in μm, some in L and some in Mbarn)
NA = 6.022e23;
κ_simple_units = 1.0e-37*NA; # simplified model
κ_units        = 1.0e-25*NA; # original model

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

#TODO: fix unit error (the value for ρ_gas_0 is good, but the units don't make sense)
M_water = 18.02e-3 # kg/mol
R_gas   = 8.314   # J/(K*mol)
T_gas   = 280.0   # K
P_gas   = 1.0e-2*1.0e5   # Pa (0.5/760.0 =>0.5 Torr?)
P_H2O   = 0.025;
ρ_gas_0 = P_H2O*P_gas*M_water/(R_gas*T_gas); # 0.1935e-3 mol m^{-3}
ρ_gas   = ρ_gas_0*(μ0./r_gas); # assume that the gas concentration profile is decreasing inversly proportional to the distance to the sample (the amount of mater crossing the surface of the cylinder is the same that crossing a bigger virtual cylinder)


Δtr = 0.5 # 1.0; # 
ρA_1 = ρ_bulk*logistic.(1000.0*(μ0.-r)/Δtr,ρ_vac,ρ0,1.0);
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

# beamprofile spread
σx = 0.5*100.0; # spread of the beam
σy = 0.5*25.0;

##
## limits of the photon energy range
##
hν1 =  897.0;
hν2 = 2117.0;

# photon flux variation w.r.t. photon energy
hknot  = [897.0; 1565.0; 1154.0; 2117.0];    
Fknot  = [96551987962033.5; 180756074879052.0; 174706655794778.0; 26025420979116.2]
μF_var = (inv([hknot.^3 hknot.^2 hknot hknot.^0])*Fknot)';

## 
## cross section density model: just a dummy model looking like O1s
## 
dKe = 0.05;
Be = collect(526.0:dKe:542.0);
Nspectrum = length(Be);

BeO1s    = 534.0; # reference binding energy of O1s under some circumstencies...
BeO1s_th = 547.0; # threshold for the background in the case of O1s (I took it away from the reference binding energy contrary to C1s)
ΔBeO1s   = 3.0;   # width of the background transition (width around the threshold)

hknot_peak = [650.0;  651.0;  789.0;      897.0; 907.0;      1154.0;   1197.0; 1565.0; 1900.0];
Eb_knot_liq = [533.84;  533.882; 534.392; 534.456; 534.581; 534.869; 534.667; 533.989; 0.998*533.989];
Eb_knot_gas = [535.818; 535.756; 536.273; 536.324; 536.438; 536.680; 536.479; 535.765; 0.998*535.765];

σEb_knot_liq = 0.5*[1.091; 1.008 ; 1.498; 1.608; 1.618; 1.624; 1.659; 1.641; 0.7*1.641];
σEb_knot_gas = 0.5*[0.629; 0.631; 0.804; 0.853;  0.863; 0.949; 0.995; 1.016; 0.7*1.016];

Ahν = [hknot_peak.^3 hknot_peak.^2 hknot_peak.^1 hknot_peak.^0];

μBe_var_liq = (inv(Ahν'*Ahν)*Ahν'*Eb_knot_liq)';
μBe_var_gas = (inv(Ahν'*Ahν)*Ahν'*Eb_knot_gas)';

σBe_var_liq = (inv(Ahν'*Ahν)*Ahν'*σEb_knot_liq)';
σBe_var_gas = (inv(Ahν'*Ahν)*Ahν'*σEb_knot_gas)';

function σ_cs_O1s_liq(hν::Cdouble,Ke::Array{Cdouble,1})
    Be_liq = hν.-Ke;
    # make the peaks vary in location and width
    μBe_O1s_liq = μBe_var_liq*[hν^3; hν^2; hν; 1.0];
    σBe_O1s_liq = σBe_var_liq*[hν^3; hν^2; hν; 1.0];
    # partial cross section (one for each chemical state)
    σ_peak = (1.0/sqrt(2.0π*σBe_O1s_liq^2))*exp.(-(Be_liq.-μBe_O1s_liq).^2/(2.0σBe_O1s_liq^2));

    # cross section value
    σ_O1s_exp(hν)*σ_peak,μBe_O1s_liq,σBe_O1s_liq
end

function σ_cs_O1s_gas(hν::Cdouble,Ke::Array{Cdouble,1})
    Be_gas = hν.-Ke;
    # make the peaks vary in location and width
    μBe_O1s_gas = μBe_var_gas*[hν^3; hν^2; hν; 1.0];
    σBe_O1s_gas = σBe_var_gas*[hν^3; hν^2; hν; 1.0];
    # partial cross section (one for each chemical state)
    σ_peak = (1.0/sqrt(2.0π*σBe_O1s_gas^2))*exp.(-(Be_gas.-μBe_O1s_gas).^2/(2.0σBe_O1s_gas^2));

    # cross section value
    σ_O1s_exp(hν)*σ_peak,μBe_O1s_gas,σBe_O1s_gas
end


# photon energy and attenuation length ranges
hν   = collect(LinRange(hν1, hν2,Ndata));  # central photon energy for each measurement
λe_s = 1.0e-3λe_exp.(hν.-BeO1s);           # attenuation length range

# for each offset between the beam center and the LJ center
for kc in [0; 1; 2; 3] # collect(0.0:0.1:3.5) #  
    # path to save location
    local save_folder = "./";
    save_folder = string(save_folder,"cylinder_radius_",μ0,"/peak_shift/")
    save_folder = string(save_folder,"eal_",Ndata,"_restricted_range/")

    ##
    ## alignment
    ##
    xc = kc*σx;
    yc = 70.0 # 75.0; # 5.0*σy; # 6.0*σy #WARNING: this value accounts for the 
    save_folder = string(save_folder,"offcenter_",kc,"/")
    # beam profile
    bp = beamProfile(xc,yc,σx,σy);

    ##
    ## acqusition parameters
    ##
    θ_aperture = 0.5*π/4                                                       # aperture's cone angle
    α_Ω        = 4π*sin(θ_aperture/2.0)^2;                                     # aperture's solid angle
    Tj         = α_Ω*ones(Cdouble,Ndata);                                      # transmission factors
    σ_ke       = 2.0*dKe*ones(Cdouble,Ndata);                                  # kinetic energy bandwidths of the analyzer (one per photon energy)
    dhν        = hν.*((1.0/25000.0)*(hν.<500.0) + (1.0/15000.0)*(hν.>=500.0)); # bandwidth of the photon beam
    Fνj        = dropdims(μF_var*[hν.^3 hν.^2 hν.^1 hν.^0]',dims=1);           # flux densities

    # dictionary where to push the data and geometry factor
    global dictAllData  = Dict()
    global dictAllGeom  = Dict()
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
            reverse(Be0),reverse(σ_cs_liq),σ_bg(μKe)*σ_bg_density(Keij,hν[j]-BeO1s_th,ΔBeO1s).+0.001);

        ##
        ## geometry factors
        ##
        local H_deom,_,H_geom,_,_,_,_,al_liq = alignmentParameter(bp,r,θ,y,x0,y0,z0,μ0,λe_s[j])
        al_liq = 1.0e-12*al_liq*κ_units/κ_simple_units # units conversion so that the unit is in μm^{-2}

        ##
        ## for the signal without noise
        ##
        local SbgO1s      = Δt*κ_units*(H_geom'*ρA_1)*Sj_bg
        local SpectrumA_1 = Δt*κ_units*(H_geom'*ρA_1)*Sj

        ##
        ## gas phase signal
        ##
        local H_deom_gas,_,H_geom_gas,_,_,_,_,al_gas = alignmentParameter(bp,r_gas,θ_gas,y_gas,x0,y0,z0,μ0,λe_s[j])
        al_gas = 1.0e-12*al_gas*κ_units/κ_simple_units # units conversion so that the unit is in μm^{-2}
        
        local Sj_gas,_,_,_ = simulateSpectrum(Fνj[j],hν[j],dhν[j],
            Keij,σ_ke[j],Tj[j],Kdisc,
            reverse(Be0),reverse(σ_cs_gas),zeros(Cdouble,Nspectrum));
        local SpectrumA_1_gas = Δt*κ_units*(H_geom_gas'*ρ_gas)*Sj_gas;

        ##
        ## add noise
        ##
        local SO1snoise = countElectrons(SbgO1s+SpectrumA_1+SpectrumA_1_gas,D_dark*Δt)


        ##
        ## push data to dicts
        ##
        local dictData = Dict( "Ke" => Keij, "Be" => Be0, "μKe" => μKe,
            "σ_cs_dens" => σ_cs_liq./σ_O1s_exp(hν[j]), "σ_cs_dens_gas" => σ_cs_gas./σ_O1s_exp(hν[j]), "σ_tot" => σ_O1s_exp(hν[j]),
            "SpectrumA_1" => SpectrumA_1, "SpectrumA_1_gas" => SpectrumA_1_gas, "Sbg" => SbgO1s, 
            "Stot" => SbgO1s+SpectrumA_1, "Snoisy" => SO1snoise,
            "T" => Tj[j], "λ" => 1.0e3λe_s[j], "F" => Fνj[j], "hν" => hν[j], "Δt" => Δt, "α" => sum(SpectrumA_1)/sum(SpectrumA_1_gas));
        local dictMetaData = Dict("μKe" => μKe,
            "σ_tot" => σ_O1s_exp(hν[j]), "peak_mode_liq" => μBe_liq, "peak_width_liq"=>σ_be_liq, "peak_mode_gas" => μBe_gas, "peak_width_gas"=>σ_be_gas,
            "T" => Tj[j], "λ" => 1.0e3λe_s[j], "F" => Fνj[j], "hν" => hν[j], "Δt" => Δt);

        local dictGeom = Dict("model" => "sharp edge cylinder + outside vapor", 
                    "hν" => hν[j], "λ" => 1.0e3λe_s[j], "radius" => μ0, "max_depth" => k0*λe0,
                        "x0" => x0, "y0" => y0, "z0" => z0, "δr" => δr, "xc" => xc, "yc" => yc, "σx" => σx, "σy" => σy, "α" => al_liq,
                        "r" => r, "r_gas" => r_gas, "H" => H_deom, "H_true" => H_geom, "H_gas"=>H_deom_gas, "H_gas_true"=>H_geom_gas, 
                        "ρ" => ρA_1, "ρ_gas" => ρ_gas)

        local dfGeom     = DataFrame(dictGeom);
        local dfData     = DataFrame(dictData);
        local dfMetaData = DataFrame(dictMetaData);

        dictAllData[Symbol(string("hν_",Int64(round(hν[j]))))]      = (eachcol(dfData),names(dfData))
        dictAllData[Symbol(string("meta_hν_",Int64(round(hν[j]))))] = (eachcol(dfMetaData),names(dfMetaData))
        dictAllGeom[Symbol(string("λe_",string(1.0e3λe_s[j])))]       = (eachcol(dfGeom),names(dfGeom))
    end

    # save_folder = string(save_folder,"new_bg/")
    save_folder = string(save_folder,"new_eal/")
    if SAVE_DATA
        mkpath(string(save_folder,exp_tag,"/"));
        XLSX.writetable(string(save_folder,exp_tag,"/data.xlsx"); dictAllData...) 
        XLSX.writetable(string(save_folder,exp_tag,"/model.xlsx"); dictAllGeom...)
    end
    global plot_sym=:Be
    include("plotData.jl")
    if SAVE_FIG
        savefig(string(save_folder,exp_tag,"/full_measurement_model.png"))
        savefig(string(save_folder,exp_tag,"/full_measurement_model.pdf"))
    end
end

