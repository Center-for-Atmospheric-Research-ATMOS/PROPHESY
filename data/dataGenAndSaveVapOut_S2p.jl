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
# using utilsFun  # for the softMax functions

# modeling XPS
using XPSpack # experiment model (geometry factor and cross section estimation)
using ATTIRE  # kinetic energy analyzer



# tags
# SHORT_RANGE = true              # select either wide range of attenuation lengths (false) or a restricted range more similar to experimental setup (true)

MODEL_5   = true               # select the number of attenuation lengths probed
MODEL_10  = false
MODEL_20  = false

save_folder = "./";
SAVE_DATA = true   # flag for saving data and model
SAVE_FIG  = true   # save the simulated data or not
SWEEPS_ON = true
nb_sweeps = 50;

# geometry setup
λe0 = 2.0e-3;        # reference penetration depth in μm
δr = 2.0e-3          # transition to vacuum layer thickness (let's set about 1 nm)
k0 = 10;             # compute the measurement model over a distance of k0*λe0
Nr = 201;            # number of discretization points in the radial dimension
Nr_lowres = 101 ;    # for the low resolution model
Nθ = 256;            # number of discretization points in the polar angle dimension
Ny = 257;            # number of discretization points in the cylinder axis dimension
μ0 = 10.0; # 20.0;   # radius of the cylinder
L = 50.0 # 60.0 # 20; # 100.0;     # height of the irradiated sample (the vertical extent of the beam is more like 20μm instead of 100μm)
x0 = sqrt(2.0)*5000.0 # (x0,y0,z0) are the coordinates of the analyzer's apperture
y0 = 0.0;
z0 = 5000.0;

# spacial discretization 
r = collect(range(μ0-k0*λe0,μ0+δr,length=Nr));
r_lowres = collect(range(μ0-k0*λe0,μ0+δr,length=Nr_lowres));
θ0 = atan(x0,z0)
θ = collect(range(θ0-π/2.0,θ0+π/2.0,Nθ));
y = collect(range(-L/2.0,L/2.0,length=Ny));

# integration inteval
Δt = 60.0 # [s]
D_dark = 0.0 # dark current coefficient is a number of electron per unit time [# s^{-1}]

# concentration profiles (4 different cases)
ρ0 = 1.0
ρ_vac = 0.0

# bulk molar concentration
# ρ_bulk =  5.0e-3; 
ρ_bulk = 10.0e-3; # 10 
Δtr = 0.5 # 1.0; # 

# unit conversion constant (some of the quantities are in μm, some in L and some in Mbarn)
NA = 6.022e23;
κ_simple_units = 1.0e-37*NA; # simplified model
κ_units        = 1.0e-25*NA; # original model

# beamprofile spread
σx = 0.5*100.0; # spread of the beam
σy = 0.5*25.0;

##
## limits of the photon energy range
##
hν1 =  532.0;
hν2 = 1726.0;

##
## variation model for binding energy, spread and photon flux w.r.t. photon energy (spline interpolation of experimental values)
##

hknot = [532.0; 789.0; 1197.0; 1726.0]
μBe1_knot = [168.996; 169.914; 169.777; 168.174]
μBe2_knot = [170.156; 171.074; 170.937; 169.334]
# σLG_knot  = 0.5e-3*[1517.535; 1471.508; 1554.653; 1340.764]; # G and L
σG_knot  = 0.5e-3*[732.133; 891.293; 632.185; 1294.616];
Fknot  = [6.5249e13; 8.7439e13; 3.5089e14; 3.476e14];

μBe1_var =  (inv([hknot.^3 hknot.^2 hknot hknot.^0])*μBe1_knot)';
μBe2_var =  (inv([hknot.^3 hknot.^2 hknot hknot.^0])*μBe2_knot)';
# μσBe_var =  (inv([hknot.^3 hknot.^2 hknot hknot.^0])*σLG_knot)';
μσBe_var =  (inv([hknot.^3 hknot.^2 hknot hknot.^0])*σG_knot)';
μF_var =  (inv([hknot.^3 hknot.^2 hknot hknot.^0])*Fknot)';

## 
## cross section density model: just a dummy model looking like S2p
## 

dKe = 0.05;                    # kinetic energy increment
Be = collect(165.0:dKe:174.5); # binding energy sampled when measuring
Nspectrum = length(Be);

BeS2p    = 169.2;  # reference binding energy of S2p under some circumstencies...
BeS2p_th = 150.0;  # threshold for the background in the case of S2p (I took it away from the reference binding energy contrary to C1s)
ΔBeS2p = 6.0;      # width of the background transition (width around the threshold)

# The width may change, let's say, from 1000 to 1300 a.u on depth range 1,5...3,7 nm, for example. Peak positions may vary up to 4-5 nm, in my case it's from 284 to 280 eV for a SDS head group on depth range mentioned above. 
function σ_cs(hν::Cdouble,Ke::Array{Cdouble,1})
    Be_cs = hν.-Ke;
    # make the peaks vary in location and width
    μBe1_val = μBe1_var*[hν^3; hν^2; hν; 1.0]
    μBe2_val = μBe2_var*[hν^3; hν^2; hν; 1.0]
    σ_be_val = μσBe_var*[hν^3; hν^2; hν; 1.0]
    # partial cross section (one for each chemical state)
    σ_peak_1 = (1.0/sqrt(2.0π*σ_be_val^2))*exp.(-(Be_cs.-μBe1_val).^2/(2.0σ_be_val^2));
    σ_peak_2 = (1.0/sqrt(2.0π*σ_be_val^2))*exp.(-(Be_cs.-μBe2_val).^2/(2.0σ_be_val^2));
    # quantity of chemical states
    p1 = 0.66 
    p2 = 0.34

    # cross section value (only for hν ∈ [295,1500.0])
    σ_S2p_exp(hν)*(p1*σ_peak_1+p2*σ_peak_2),[μBe1_val;μBe2_val],σ_be_val,[p1;p2]
end

FLAG_OFF_CENTER = [true false false false;
                    false true false false;
                    false false true false;
                    false false false true];
FLAG_PROFILE    = [true false false false false;
                    false true false false false;
                    false false true false false;
                    false false false true false;
                    false false false false true];

# FLAG_PROFILE = [false false false false true]';

# for (FLAG_OFF_CENTER_0,FLAG_OFF_CENTER_1,FLAG_OFF_CENTER_2,FLAG_OFF_CENTER_3,FLAG_0001,FLAG_0002,FLAG_0003,FLAG_0004) in zip(FLAG_OFF_CENTER,FLAG_PROFILE)
for (FLAG_OFF_CENTER_0,FLAG_OFF_CENTER_1,FLAG_OFF_CENTER_2,FLAG_OFF_CENTER_3) in eachcol(FLAG_OFF_CENTER)
    for (FLAG_0001,FLAG_0002,FLAG_0003,FLAG_0004,FLAG_0005) in eachcol(FLAG_PROFILE)
        # soon add the option of sphere or plane geometry?
        local save_folder = "./";
        save_folder = string(save_folder,"cylinder_radius_",μ0,"/peak_shift/")

        if FLAG_0001
            ρA_1 = ρ_bulk*logistic.(1000.0*(μ0.-r)/Δtr,ρ_vac,ρ0,1.0);
            global exp_tag     = "0001/S2p"
        end
        if FLAG_0002
            ρA_1 = ρ_bulk*(logistic.(1000.0*(μ0.-r)/Δtr,ρ_vac,ρ0,1.0) .+ 2.0exp.(-(1000.0*(μ0.-r).-0.0).^2. /(2.0*(0.5Δtr)^2)));
            global exp_tag     = "0002/S2p"
        end
        if FLAG_0003
            ρA_1 = ρ_bulk*(logistic.(1000.0*(μ0.-r)/Δtr,ρ_vac,ρ0,1.0) .+ exp.(-(1000.0*(μ0.-r).-0.5).^2. /(2.0*Δtr^2)));
            global exp_tag     = "0003/S2p"
        end
        if FLAG_0004
            ρA_1 = ρ_bulk*exp.(-(1000.0((μ0.-r))).^2. /(2.0*Δtr^2));
            global exp_tag     = "0004/S2p"
        end
        if FLAG_0005
            ρA_1 = ρ_bulk*(logistic.(1000.0*(μ0.-r)/Δtr,ρ_vac,ρ0,1.0) .+ exp.(-(1000.0*(μ0.-r).-0.5).^2. /(2.0*(2.0Δtr)^2)));
            global exp_tag     = "0005/S2p"
        end

        # measurement operator (only the geometrical term since the other comes as a multiplicative scalar estimated from the data)
        global  Ndata = 5;
        if MODEL_5                                   # number of measurement (penetration depth)
            Ndata = 5;
        end
        if MODEL_10
            Ndata = 10;
        end
        if MODEL_20
            Ndata = 20;
        end


        save_folder = string(save_folder,"eal_",Ndata,"_restricted_range/")
        global hν = collect(LinRange(hν1, hν2,Ndata));                # central photon energy for each measurement
        global λe = 1.0e-3λe_exp.(hν.-BeS2p);              # attenuation length range

        ##
        ## alignment
        ##
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
        xc = kc*σx; # alignment with the photon flux
        yc = 70.0 # 75.0; # 99.0; # 5.0*σy; # 6.0*σy #WARNING: this value accounts for the misalignment with the analyzer

        save_folder = string(save_folder,"offcenter_",kc,"/")
        # beam profile
        bp = beamProfile(xc,yc,σx,σy);

        
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
        Fνj        = dropdims(μF_var*[hν.^3 hν.^2 hν.^1 hν.^0]',dims=1);

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
            local σ_cs_fg,μBe_dev,σ_be_var,τ_be =  σ_cs(hν[j],Keij);

            ##
            ## electron flux without the geometry factor: signal of interest (Sj) and background signal (Sj_bg)
            ##
            # local Sj,Sj_bg,_,_ = simulateSpectrum(Fνj[j],hν[j],dhν[j],
            #     Keij,σ_ke[j],Tj[j],Kdisc,
            #     reverse(Be0),reverse(σ_cs_fg),σ_bg(μKe)*σ_bg_density(Keij,hν[j]-BeS2p_th,ΔBeS2p)); # σ_bg_lin_density(Keij,hν[j]-BeS2p_th,5.0e4ΔBeC1s)
            local Sj,Sj_bg,_,_ = simulateSpectrum(Fνj[j],hν[j],dhν[j],
                Keij,σ_ke[j],Tj[j],Kdisc,
                reverse(Be0),reverse(σ_cs_fg),25.0*σ_bg(μKe)*σ_bg_density(Keij,hν[j]-BeS2p_th,ΔBeS2p)); # σ_bg_lin_density(Keij,hν[j]-BeS2p_th,5.0e4ΔBeC1s) 5.0

            ##
            ## geometry factors
            ##
            local H_deom,_,H_geom,_,_,_,_,al_C1s = alignmentParameter(bp,r,θ,y,x0,y0,z0,μ0,λe[j])
            al_C1s = 1.0e-12*al_C1s*κ_units/κ_simple_units # units conversion so that the unit is in μm^{-2}

            ##
            ## the signal without noise (time integrated)
            ##
            local SbgC1s      = Δt*κ_units*(H_geom'*ρA_1)*Sj_bg
            local SpectrumA_1 = Δt*κ_units*(H_geom'*ρA_1)*Sj

            ##
            ## add noise
            ##
            
            local SC1snoise = countElectrons(SbgC1s+SpectrumA_1,D_dark*Δt)

            if SWEEPS_ON # just repeat the acquisition
                local SC1snoise_sweeps = zeros(Cdouble,length(SC1snoise),nb_sweeps)
                for i in 1:nb_sweeps
                    SC1snoise_sweeps[:,i] = countElectrons(SbgC1s+SpectrumA_1,D_dark*Δt);
                end
            end

            # plot signals w.r.t. the kinetic energy
            # figure(); plot(Keij,SbgC1s); plot(Keij,SbgC1s+SpectrumA_1); scatter(Keij,rand.(Poisson.(SbgC1s+SpectrumA_1))); xlabel("kinetic energy [eV]"); ylabel("spectrum [a.u.]") 

            ##
            ## push data to dicts
            ##

            local dictData = Dict( "Ke" => Keij, "Be" => Be0, "μKe" => μKe,
                "σ_cs_dens" => σ_cs_fg./σ_S2p_exp(hν[j]), "σ_tot" => σ_S2p_exp(hν[j]),
                "SpectrumA_1" => SpectrumA_1, "Sbg" => SbgC1s, 
                "Stot" => SbgC1s+SpectrumA_1, "Snoisy" => SC1snoise,
                "T" => Tj[j], "λ" => 1.0e3λe[j], "F" => Fνj[j], "hν" => hν[j], "Δt" => Δt);
            if SWEEPS_ON
                [dictData[string("Snoisy_",i)] = SC1snoise_sweeps[:,i] for i in 1:nb_sweeps]
            end
            local dictMetaData = Dict("μKe" => μKe,
                "σ_tot" => σ_S2p_exp(hν[j]), "peak_mode" => μBe_dev, "peak_width"=>σ_be_var, "peak_probability"=>τ_be,
                "T" => Tj[j], "λ" => 1.0e3λe[j], "F" => Fνj[j], "hν" => hν[j], "Δt" => Δt);

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

        # save_folder = string(save_folder,"new_bg/")
        save_folder = string(save_folder,"new_eal/")
        if SAVE_DATA
            mkpath(string(save_folder,exp_tag,"/"));
            XLSX.writetable(string(save_folder,exp_tag,"/data.xlsx"); dictAllData...) # TODO: get outside the loop
            XLSX.writetable(string(save_folder,exp_tag,"/model.xlsx"); dictAllGeom...)
        end
        global plot_sym=:Be
        include("plotData.jl")
        if SAVE_FIG
            savefig(string(save_folder,exp_tag,"/full_measurement_model.png"))
            savefig(string(save_folder,exp_tag,"/full_measurement_model.pdf"))
        end
    end
end

