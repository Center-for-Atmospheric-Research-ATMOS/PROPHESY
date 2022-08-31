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
# SHORT_RANGE = true              # select either wide range of attenuation lengths (false) or a restricted range more similar to experimental setup (true)

MODEL_5   = true               # select the number of attenuation lengths probed
MODEL_10  = false
MODEL_20  = false

# FLAG_OFF_CENTER_0 = true;
# FLAG_OFF_CENTER_1 = false;
# FLAG_OFF_CENTER_2 = false;
# FLAG_OFF_CENTER_3 = false;

# FLAG_0001 = false               # selection of the profile (one must be true and the others false)
# FLAG_0002 = false
# FLAG_0003 = true
# FLAG_0004 = false

save_folder = "./";
SAVE_DATA = false   # flag for saving data and model
SAVE_FIG  = false   # save the simulated data or not
SWEEPS_ON = false
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

# integration inteval
Δt = 60.0 # [s]
D_dark = 0.0 # dark current coefficient is a number of electron per unit time [# s^{-1}]

# unit conversion constant (some of the quantities are in μm, some in L and some in Mbarn)
NA = 6.022e23;
κ_simple_units = 1.0e-37*NA; # simplified model
κ_units        = 1.0e-25*NA; # original model

FLAG_OFF_CENTER = [true false false false;
                    false true false false;
                    false false true false;
                    false false false true];
FLAG_PROFILE    = [true false false false;
                    false true false false;
                    false false true false;
                    false false false true];

# for (FLAG_OFF_CENTER_0,FLAG_OFF_CENTER_1,FLAG_OFF_CENTER_2,FLAG_OFF_CENTER_3,FLAG_0001,FLAG_0002,FLAG_0003,FLAG_0004) in zip(FLAG_OFF_CENTER,FLAG_PROFILE)
for (FLAG_OFF_CENTER_0,FLAG_OFF_CENTER_1,FLAG_OFF_CENTER_2,FLAG_OFF_CENTER_3) in eachcol(FLAG_OFF_CENTER)
    for (FLAG_0001,FLAG_0002,FLAG_0003,FLAG_0004) in eachcol(FLAG_PROFILE)
        # soon add the option of sphere or plane geometry?
        local save_folder = "./";
        save_folder = string(save_folder,"cylinder_radius_",μ0,"/peak_shift/")

        # spacial discretization 
        r = collect(range(μ0-k0*λe0,μ0+δr,length=Nr));
        r_lowres = collect(range(μ0-k0*λe0,μ0+δr,length=Nr_lowres));
        θ0 = atan(x0,z0)
        θ = collect(range(θ0-π/2.0,θ0+π/2.0,Nθ));
        y = collect(range(-L/2.0,L/2.0,length=Ny));

        # concentration profiles (4 different cases)
        ρ0 = 1.0
        ρ_vac = 0.0

        # bulk molar concentration
        # ρ_bulk =  5.0e-3; 
        ρ_bulk = 10.0e-3; # 10 
        Δtr = 0.5 # 1.0; # 


        if FLAG_0001
            ρA_1 = ρ_bulk*logistic.(1000.0*(μ0.-r)/Δtr,ρ_vac,ρ0,1.0);
            global exp_tag     = "0001"
        end
        if FLAG_0002
            ρA_1 = ρ_bulk*(logistic.(1000.0*(μ0.-r)/Δtr,ρ_vac,ρ0,1.0) .+ 2.0exp.(-(1000.0*(μ0.-r).-0.0).^2. /(2.0*(0.5Δtr)^2)));
            global exp_tag     = "0002"
        end
        if FLAG_0003
            ρA_1 = ρ_bulk*(logistic.(1000.0*(μ0.-r)/Δtr,ρ_vac,ρ0,1.0) .+ exp.(-(1000.0*(μ0.-r).-0.5).^2. /(2.0*Δtr^2)));
            global exp_tag     = "0003"
        end
        if FLAG_0004
            ρA_1 = ρ_bulk*exp.(-(1000.0((μ0.-r))).^2. /(2.0*Δtr^2));
            global exp_tag     = "0004"
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

        # not sure it's that short anymore
        # if SHORT_RANGE                                          # bounds of the attenuation lengths
            # λe1 = 1.3; hν1 = 365.0;
            # λe2 = 2.5; hν2 = 1500.0;
            λe1 = 1.3; hν1 =  650.0;
            λe2 = 3.8; hν2 = 1884.0;
            save_folder = string(save_folder,"eal_",Ndata,"_restricted_range/")
        # else
        #     λe1 = 0.5; hν1 = 310.0;
        #     λe2 = 5.5; hν2 = 1900.0;
        #     save_folder = string(save_folder,"eal_",Ndata,"/")
        # end
        # global λe = 1.0e-3collect(range(λe1,λe2,Ndata));              # attenuation length range
        global hν = collect(LinRange(hν1, hν2,Ndata));                # central photon energy for each measurement
        BeC1s = 284.0;
        global λe = 1.0e-3λe_exp.(hν.-BeC1s);              # attenuation length range

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
        yc = 70.0 # 75.0; # 5.0*σy; # 6.0*σy #WARNING: this value accounts for the 

        save_folder = string(save_folder,"offcenter_",kc,"/")
        # beam profile
        bp = beamProfile(xc,yc,σx,σy);
        # α_al = zeros(Cdouble,Ndata);
        # for i in 1:Ndata
        #     local H_r,H_rθy,H_r_ph,H_rθy_ph,_,_,_,α_al[i] =  alignmentParameter(bp,r,θ,y,x0,y0,z0,μ0,λe[i]);
        # end

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

        BeC1s_th = 284.0 # mean(Be);
        ΔBeC1s = 6.0 # 12.0/5.0#(Be[end]-Be[1])/5;

        # The width may change, let's say, from 1000 to 1300 a.u on depth range 1,5...3,7 nm, for example. Peak positions may vary up to 4-5 nm, in my case it's from 284 to 280 eV for a SDS head group on depth range mentioned above. 

        Ahν  = [650.0^2  650.0 1.0; 
        1315.0^2 1315.0 1.0;
        1884.0^2 1884.0 1.0];
        Bhν = [1.0; (286.0+1.0)/286.0; (286.0-3.0)/286.0];
        μBe_var = inv(Ahν)*Bhν;

        hhν = collect(310.0:20.0:1900.0);
        # figure()
        # scatter([650.0; 1315.0; 1884.0],Bhν)
        # plot(hhν,dropdims(μBe_var'*[hhν.^2 hhν.^1 hhν.^0]',dims=1))


        # figure();
        # plot(collect(650.0:10.0:1884.0),dropdims(μBe_var'*[collect(650.0:10.0:1884.0).^2 collect(650.0:10.0:1884.0) ones(Cdouble,length(collect(650.0:10.0:1884.0)))]',dims=1))

        function σ_cs(hν::Cdouble,Ke::Array{Cdouble,1},μKe::Cdouble;μKe0::Cdouble=50.0,μKe1::Cdouble=1200.0)
            Be_cs = hν.-Ke;
            # make the peaks vary in location and width
            peak_dev = μBe_var'*[hν^2; hν; 1.0];
            μBe_dev = peak_dev*μBe
            # σ_be_var = (1.0+0.3*((μKe-μKe0)/(μKe1-μKe0)))*σ_be;
            σ_be_var = (1.0+0.3*((hν-650.0)/(1884.0-650.0)))*σ_be;
            # partial cross section (one for each chemical state)
            σ_peak_1 = (1.0/sqrt(2.0π*σ_be_var[1]^2))*exp.(-(Be_cs.-μBe_dev[1]).^2/(2.0σ_be_var[1]^2));
            σ_peak_2 = (1.0/sqrt(2.0π*σ_be_var[2]^2))*exp.(-(Be_cs.-μBe_dev[2]).^2/(2.0σ_be_var[2]^2));
            σ_peak_3 = (1.0/sqrt(2.0π*σ_be_var[3]^2))*exp.(-(Be_cs.-μBe_dev[3]).^2/(2.0σ_be_var[3]^2));
            # quantity of chemical states
            p1 = 0.7 # 0.85 .+ (0.77-0.85)*(μKe.-μKe0)./(μKe1-μKe0); # 0.8 # 
            p2 = 0.25 # 0.125 .+ (0.12-0.125)*(μKe.-μKe0)./(μKe1-μKe0); # 0.12 # 
            p3 = 1.0-(p1+p2);

            # cross section value (only for hν ∈ [295,1500.0])
            σ_C1s_exp(hν)*(p1*σ_peak_1+p2*σ_peak_2+p3*σ_peak_3),μBe_dev,σ_be_var,[p1;p2;p3]
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
            # local σ_cs_fg,μBe_dev,σ_be_var,τ_be =  σ_cs.(hν[j],Keij,μKe);
            local σ_cs_fg,μBe_dev,σ_be_var,τ_be =  σ_cs(hν[j],Keij,μKe);

            ##
            ## electron flux without the geometry factor: signal of interest (Sj) and background signal (Sj_bg)
            ##
            # local Sj,Sj_bg,_,_ = simulateSpectrum(Fνj[j],hν[j],dhν[j],
            #     Keij,σ_ke[j],Tj[j],Kdisc,
            #     reverse(Be0),reverse(σ_cs_fg),σ_bg(μKe)*σ_bg_density(Keij,hν[j]-BeC1s_th,ΔBeC1s)); # σ_bg_lin_density(Keij,hν[j]-BeC1s_th,5.0e4ΔBeC1s)
            local Sj,Sj_bg,_,_ = simulateSpectrum(Fνj[j],hν[j],dhν[j],
                Keij,σ_ke[j],Tj[j],Kdisc,
                reverse(Be0),reverse(σ_cs_fg),σ_bg(μKe)*σ_bg_density(Keij,hν[j]-BeC1s_th,ΔBeC1s)); # σ_bg_lin_density(Keij,hν[j]-BeC1s_th,5.0e4ΔBeC1s) 5.0

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
                "σ_cs_dens" => σ_cs_fg./σ_C1s_exp(hν[j]), "σ_tot" => σ_C1s_exp(hν[j]),
                "SpectrumA_1" => SpectrumA_1, "Sbg" => SbgC1s, 
                "Stot" => SbgC1s+SpectrumA_1, "Snoisy" => SC1snoise,
                "T" => Tj[j], "λ" => 1.0e3λe[j], "F" => Fνj[j], "hν" => hν[j], "Δt" => Δt);
            if SWEEPS_ON
                [dictData[string("Snoisy_",i)] = SC1snoise_sweeps[:,i] for i in 1:nb_sweeps]
            end
            local dictMetaData = Dict("μKe" => μKe,
                "σ_tot" => σ_C1s_exp(hν[j]), "peak_mode" => μBe_dev, "peak_width"=>σ_be_var, "peak_probability"=>τ_be,
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


if false

    # [a] E. Kukk, J.D. Bozek, G. Snell, W.-T. Cheng and N. Berrah, Phys. Rev. A v.63, 062702 (2001).
    # [b] E. Kukk, K. Ueda, U. Hergenhahn, J. Liu X, G. Prumper, H. Yoshida, Y. Tamenori, C. Makochekanwa, T. Tanaka, M. Kitajima and H. Tanaka, Phys.Rev.Lett. v. 95, p. 133001 (2005).

    τm = [0.85; 0.125; 1.0-0.85-0.125];  # [1.0/3.0; 1.0/3.0; 1.0/3.0];
    μm = [290.2; 292.0; 293.0]; # [290.3; 291.9; 293.5];
    μm = [280.0; 281.0; 284.0]
    σm = sqrt(2.0)*[0.45; 0.25; 0.6]; # [290.3; 291.9; 293.5]/500.0;
    τt = zeros(Cdouble,Ndata,3);
    μt = zeros(Cdouble,Ndata,3);
    σt = zeros(Cdouble,Ndata,3);
    for j in 1:Ndata
        local symbol_h = Symbol(string("hν_",Int64(round(hν[j]))));
        local be = dictAllData[symbol_h][1][:Be]; # reverse();
        # local spectrum = dictAllData[symbol_h][1][:SpectrumA_1];
        local peak_dev = μBe_var'*[hν[j]^2; hν[j]; 1.0];
        local spectrum = dictAllData[symbol_h][1][:Snoisy]-dictAllData[symbol_h][1][:Sbg];
    # estimate the peaks centers and spreads
    τt[j,:],μt[j,:],σt[j,:] = EM_peaks(be,spectrum,τm,peak_dev*μm,σm,200)
    end

    figure()
    for j in  1:Ndata # 1:5:Ndata
        # local μKe0=50.0
        # local μKe1=1200.0
        local symbol_h = Symbol(string("hν_",Int64(round(hν[j]))));
        local be = dictAllData[symbol_h][1][:Be];
        local μKe = dictAllData[symbol_h][1][:μKe];
        local Keij = reverse(hν[j] .- Be) ; 
    
        local σ_peak,_,_,_ =  σ_cs(hν[j],Keij,μKe[1])
        σ_peak = σ_peak/σ_C1s_exp(hν[j])

        # estimation
        σ_est_1 = τt[j,1]*(1.0/sqrt(2.0π*σt[j,1]^2))*exp.(-(be.-μt[j,1]).^2/(2.0σt[j,1]^2));
        σ_est_2 = τt[j,2]*(1.0/sqrt(2.0π*σt[j,2]^2))*exp.(-(be.-μt[j,2]).^2/(2.0σt[j,2]^2));
        σ_est_3 = τt[j,3]*(1.0/sqrt(2.0π*σt[j,3]^2))*exp.(-(be.-μt[j,3]).^2/(2.0σt[j,3]^2));

        # println(dKe*sum(p1*σ_peak_1+p2*σ_peak_2+p3*σ_peak_3))
        println(dKe*sum(σ_peak))
        println(dKe*sum(σ_est_1+σ_est_2+σ_est_3),"\n")

        # plot(be,p1*σ_peak_1+p2*σ_peak_2+p3*σ_peak_3)
        plot(be,σ_peak)
        scatter(be,σ_est_1+σ_est_2+σ_est_3)

        # plot(be,(dictAllData[symbol_h][1][:Snoisy]-dictAllData[symbol_h][1][:Sbg])/(dKe*sum(dictAllData[symbol_h][1][:Snoisy]-dictAllData[symbol_h][1][:Sbg])))
    end
end

if false
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
        local H_r,H_rθy,H_r_ph,H_rθy_ph,_,_,_,α_al_off[ic] =  alignmentParameter(bp,r,θ,y,x0,y0,z0,μ0,λe[1]);
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
        α_al_noise[i] = τ_al_noise[i]/(Tj[i]*Fνj[i]*σ_C1s_exp(hν[i]))
        α_al_noise_gt[i] = τ_al_noise_gt[i]/(Tj[i]*Fνj[i]*σ_C1s_exp(hν[i]))
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