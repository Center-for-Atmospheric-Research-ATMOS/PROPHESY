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
FLAG_0002 = false
FLAG_0003 = false
FLAG_0004 = true

save_folder = "./";

# geometry setup
λe0 = 2.0e-3;        # reference penetration depth in μm
δr = 2.0e-3          # transition to vacuum layer thickness (let's set about 1 nm)
k0 = 10;             # compute the measurement model over a distance of k0*λe0
Nr = 201;            # number of discretization points in the radial dimension
Nr_lowres = 101 ;    # for the low resolution model
Nθ = 256;            # number of discretization points in the polar angle dimension
Ny = 256;            # number of discretization points in the cylinder axis dimension
μ0 = 10.0; # 20.0;   # radius of the cylinder
L = 100.0;           # height of the irradiated sample
x0 = sqrt(2.0)*100.0 # (x0,y0,z0) are the coordinates of the analyzer's apperture
y0 = 0.0;
z0 = 100.0;

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
r_th  = 2.346; # 2.0


if FLAG_0001
    ρA_1 = logistic.(1000.0*(δr.+μ0.-r).-r_th,ρ_vac,ρ0,2.0);
    exp_tag     = "0001"
end
if FLAG_0002
    ρA_1 = logistic.(1000.0*(δr.+μ0.-r).-r_th,ρ_vac,ρ0,2.0) .+ 2.0exp.(-(1000.0*(μ0.-r).-0.0).^2. /(2.0*0.25^2));
    exp_tag     = "0002"
end
if FLAG_0003
    ρA_1 = logistic.(1000.0*(δr.+μ0.-r).-r_th,ρ_vac,ρ0,2.0) .+ exp.(-(1000.0*(μ0.-r).-0.5).^2. /(2.0*0.5^2));
    exp_tag     = "0003"
end
if FLAG_0004
    ρA_1 = exp.(-(1000.0((δr.+μ0.-r).-δr)).^2. /(2.0*0.5^2));
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
    λe1 = 1.3
    λe2 = 2.5
    save_folder = string(save_folder,"eal_",Ndata,"_restricted_range/")
else
    λe1 = 0.5;
    λe2 = 5.5; 
    save_folder = string(save_folder,"eal_",Ndata,"/")
end
λe = 1.0e-3collect(range(λe1,λe2,Ndata));              # attenuation length range

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
## simple analyzer model
##
"""
    (Fνj::Cdouble,hνj::Cdouble,Δνj::Cdouble): photon beam's parameters
    (Ki::Array{Cdouble,1},ΔKi::Cdouble,T::Cdouble): analyzer's parameters
    (r::Array{Cdouble,1},θ::Array{Cdouble,1},y::Array{Cdouble,1},x0,y0,z0,μ0,λe): geometry factor (parameters)
    (hν::Array{Cdouble,1},Ki::Array{Cdouble,1},Be0::Array{Cdouble,1},σ_cs_0::Array{Cdouble,1}): cross section
    (ρ::Array{Cdouble,1}): concentration profile
"""
function simulateSpectrumGeom(Fνj::Cdouble,hνj::Cdouble,Δνj::Cdouble,
    Ki::Array{Cdouble,1},ΔKi::Cdouble,T::Cdouble,
    Be0::Array{Cdouble,1},σ_cs_0::Array{Cdouble,1},
    r::Array{Cdouble,1},θ::Array{Cdouble,1},y::Array{Cdouble,1},x0::Cdouble,y0::Cdouble,z0::Cdouble,μ0::Cdouble,λ::Cdouble,
    ρ::Array{Cdouble,1})

    ##
    ## analyzer
    ##
    # efficiency functions of the analyzer
    Nchannel = length(Ki); # number of readings in a spectrum
    dKe = Ki[2]-Ki[1];
    φi = Φi(Ki,ΔKi,T); 

    ##
    ## photon beam spectrum
    ##
    # photon energy discretization space
    nν = 5 # this should depend on the relative value ΔKi and Δνj
    Δhν = dKe; # discretization step in the photon energy space, note: not the same as the bandwith of the photon spectrum Δνj
    hνjl = collect(hνj-nν*Δhν:Δhν:hνj+nν*Δhν); # discretization of the photon energy space

    ##
    ## spread of the ionization cross section: light source and analyzer
    ##
    # number of discretization point in the kinetic energy space and in the photon energy space
    Nspectrum = length(Ki); # NOTE: the discretization of the kinetic enegy space does not have to coincide with the center of the channel, it can be whatever subdivision of that space
    NdensF = length(hνjl);
    # discretization of the integral: piecewise linear basis function
    Gm = dKe*ones(Cdouble,Nspectrum); Gm[1] = 0.5dKe; Gm[end] = 0.5dKe;
    Fl = Δhν*ones(Cdouble,NdensF); Fl[1] = 0.5Δhν; Fl[end] = 0.5Δhν;

    # define an interpolating tool whose nodes are Be0::Array{Cdouble,1},σ_cs_0::Array{Cdouble,1}
    σ_cs_interp = extrapolate(interpolate((Be0,), σ_cs_0, Gridded(Linear())),Line())

    # compute the "convolution" of the corss section by the spread of the light source and the spread of the kinetic energy analyzer
    Aij = zeros(Cdouble,Nspectrum,NdensF);  # discretization of the spread of the source and analyzer
    Sj = Array{Cdouble,1}(undef,Nspectrum); 
    
    F_dens = sourceSpread.(hνjl,hνj,Δνj,Fνj)
    σ_tot = XPSpack.σ_C1s_interp[hνjl] 
    σ_val = Array{Cdouble,2}(undef,NdensF,Nspectrum)
    for l in 1:NdensF
        # interpolator for the current photon energy
        Be = hνjl[l] .- Ki;
        σ_val[l,:] = σ_tot[l]*σ_cs_interp[Be] # no need to repeat that computation at each iteration of the i loop, once at first is enough
    end
    σ_val[σ_val.<0.0] .= 0.0 ;
    for i in 1:Nchannel 
        φi_val = φi[i].(Ki)
        for m in 1:Nspectrum
            for l in 1:NdensF
                Aij[m,l] = φi_val[m]*F_dens[l]*σ_val[l,m] # σ_cs(hνjl[l],Keij[m],μKe);
            end
        end
        Sj[i] = Gm'*Aij*Fl;
    end


    ##
    ## geometry factor
    ##
    H_geom,_,_,_,_ = cylinder_gain_H(r,θ,y,x0,y0,z0,μ0,λ);
    
    ##
    ## compute total signal
    ##
    (H_geom'*ρ)*Sj,H_geom,Sj,Keij
end

# acqusition parameters
# j = Ndata; # 1; # 10;                                               # select the photon energy
hν = collect(LinRange(365.0,1500.0,Ndata));                         # central photon energy for each measurement
dhν = hν.*((1.0/25000.0)*(hν.<500.0) + (1.0/15000.0)*(hν.>=500.0)); # bandwidth of the photon beam
Fνj = 1.0e3*ones(Cdouble,Ndata);                                    # flux densities
Tj   = collect(LinRange(5.0,10.0,Ndata));                           # transmission factors
σ_ke = 2.0*dKe*collect(LinRange(1.0,2.0,Ndata));                    # kinetic energy bandwidths of the analyzer (one per photon energy)

# dictionary where to push the data and geometry factor
dictAllData = Dict() # TODO: get outside the loop
dictAllGeom = Dict()
for j in 1:Ndata # can potentially multi-thread this loop: but need to sync before writing files
    # for a given photon energy νj, measure a spectrum
    local Keij = reverse(hν[j] .- Be) ;                                       # centers of the analyzer's channels
    local μKe = 0.5*(Keij[1]+Keij[end]);                                      # central kinetic energy (a bit the same role as pass energy)

    ##
    ## compute the cross section of the sample (to be estimated from the data in an estimation setting)
    ##

    local Be0 = hν[j] .- Keij;
    local σ_cs_0 =  σ_cs.(hν[j],Keij,μKe)./XPSpack.σ_C1s_interp[hν[j]] 
    local SpectrumA_1,H_geom,S_anph,Ki = simulateSpectrumGeom(Fνj[j],hν[j],dhν[j],
        Keij,σ_ke[j],Tj[j],
        reverse(Be0),reverse(σ_cs_0),
        r,θ,y,x0,y0,z0,μ0,λe[j],
        ρA_1)

    ##
    ## add background and noise
    ##

    # parameters used for the simulation of the inelastic background 
    # Here, the background is made up of electrons that undergo inelastic collisions,
    # at least one so that they don't appear in the sharp peak, but at lower kinetic energy
    # the model is fairly simple and arbitrary, but it's better than no the background at all
    local SbgC1s = 0.5Keij./(1.0.+exp.((Keij.-(hν[j]-BeC1s))./ΔBeC1s)); # the inelastic collision background

    # add the noise on the signal that hit the sensor, i.e. signal with background
    # SC1snoise = countElectrons(SbgC1s+Ssignal[:,j])
    local SC1snoise = countElectrons(SbgC1s+SpectrumA_1)

    # SC1s = SC1snoise - SbgC1s; # noisy signal without background -> negative values appears

    # plot signals w.r.t. the kinetic energy
    figure(); plot(Keij,SbgC1s); plot(Keij,SbgC1s+SpectrumA_1); scatter(Keij,rand.(Poisson.(SbgC1s+SpectrumA_1))); xlabel("kinetic energy [eV]"); ylabel("spectrum [a.u.]") 
    # plot the signal w.r.t. the binding energy
    # figure(); plot(Be,reverse(SpectrumA_1)); scatter(Be,reverse(SC1snoise-SbgC1s)); ax = gca(); ax.invert_xaxis(); xlabel("binding energy [eV]"); ylabel("spectrum (no background) [a.u.]") 


    ##
    ## push data to dicts
    ##

    local dictData = Dict( "Ke" => Keij, "Be" => Be0, "μKe" => μKe,
        "σ_cs_dens" => σ_cs_0, "σ_tot" => XPSpack.σ_C1s_interp[hν[j]], 
        "SpectrumA_1" => SpectrumA_1, "Sbg" => SbgC1s, 
        "Stot" => SbgC1s+SpectrumA_1, "Snoisy" => SC1snoise,
        "T" => Tj[j], "λ" => 1.0e3λe[j], "F" => Fνj[j], "hν" => hν[j]);

    local dictGeom = Dict("model" => "sharp edge cylinder + outside vapor", 
                "hν" => hν[j], "λ" => 1.0e3λe[j], "radius" => μ0, "max_depth" => k0*λe0,
                    "x0" => x0, "y0" => y0, "z0" => z0, "δr" => δr,
                    "r" => r, "H" => H_geom, "ρ" => ρA_1)

    local dfGeom = DataFrame(dictGeom);
    local dfData = DataFrame(dictData);

    dictAllData[Symbol(string("hν_",Int64(round(hν[j]))))] = (eachcol(dfData),names(dfData))
    dictAllGeom[Symbol(string("λe_",string(1.0e3λe[j])))] = (eachcol(dfGeom),names(dfGeom))
end

# mkpath(string(save_folder,exp_tag,"/"));
# XLSX.writetable(string(save_folder,exp_tag,"/data.xlsx"); dictAllData...) # TODO: get outside the loop
# XLSX.writetable(string(save_folder,exp_tag,"/model.xlsx"); dictAllGeom...)

figure(figsize=[12, 10]); 
ax1 = subplot(221)
title("Eph = 365 [eV]",fontsize=14)
symbol_h = :hν_365
plot(dictAllData[symbol_h][1].Ke,dictAllData[symbol_h][1].Sbg,label="background"); 
plot(dictAllData[symbol_h][1].Ke,dictAllData[symbol_h][1].Sbg+dictAllData[symbol_h][1].SpectrumA_1,label="noise free spectrum"); 
scatter(dictAllData[symbol_h][1].Ke,dictAllData[symbol_h][1].Snoisy,label="noisy spectrum"); 
# GeomGain = dictAllGeom[Symbol("λe_1.3000000000000003")][1].H'*dictAllGeom[Symbol("λe_1.3000000000000003")][1].ρ
# plot(dictAllData[symbol_h][1].Ke,dictAllData[symbol_h][1].Sbg+GeomGain*dictAllData[symbol_h][1].σ_cs_dens.*dictAllData[symbol_h][1].σ_tot.*dictAllData[symbol_h][1].T.*dictAllData[symbol_h][1].F,label="unaltered spectrum"); 
# scatter(Keij,rand.(Poisson.(SbgC1s+SpectrumA_1))); 
xlabel("kinetic energy [eV]",fontsize=14); ylabel("spectrum [a.u.]",fontsize=14) 
xticks(fontsize=14); yticks(fontsize=14); 
legend(fontsize=14)

ax2 = subplot(222)
title("Eph = 649 [eV]",fontsize=14)
symbol_h = :hν_649
plot(dictAllData[symbol_h][1].Ke,dictAllData[symbol_h][1].Sbg,label="background"); 
plot(dictAllData[symbol_h][1].Ke,dictAllData[symbol_h][1].Sbg+dictAllData[symbol_h][1].SpectrumA_1,label="noise free spectrum"); 
scatter(dictAllData[symbol_h][1].Ke,dictAllData[symbol_h][1].Snoisy,label="noisy spectrum"); 
# GeomGain = dictAllGeom[Symbol("λe_1.3000000000000003")][1].H'*dictAllGeom[Symbol("λe_1.3000000000000003")][1].ρ
# plot(dictAllData[symbol_h][1].Ke,dictAllData[symbol_h][1].Sbg+GeomGain*dictAllData[symbol_h][1].σ_cs_dens.*dictAllData[symbol_h][1].σ_tot.*dictAllData[symbol_h][1].T.*dictAllData[symbol_h][1].F,label="unaltered spectrum"); 
# scatter(Keij,rand.(Poisson.(SbgC1s+SpectrumA_1))); 
xlabel("kinetic energy [eV]",fontsize=14); ylabel("spectrum [a.u.]",fontsize=14) 
xticks(fontsize=14); yticks(fontsize=14); 
legend(fontsize=14)

ax3 = subplot(223)
title("Eph = 932 [eV]",fontsize=14)
symbol_h = :hν_932
plot(dictAllData[symbol_h][1].Ke,dictAllData[symbol_h][1].Sbg,label="background"); 
plot(dictAllData[symbol_h][1].Ke,dictAllData[symbol_h][1].Sbg+dictAllData[symbol_h][1].SpectrumA_1,label="noise free spectrum"); 
scatter(dictAllData[symbol_h][1].Ke,dictAllData[symbol_h][1].Snoisy,label="noisy spectrum"); 
# GeomGain = dictAllGeom[Symbol("λe_1.3000000000000003")][1].H'*dictAllGeom[Symbol("λe_1.3000000000000003")][1].ρ
# plot(dictAllData[symbol_h][1].Ke,dictAllData[symbol_h][1].Sbg+GeomGain*dictAllData[symbol_h][1].σ_cs_dens.*dictAllData[symbol_h][1].σ_tot.*dictAllData[symbol_h][1].T.*dictAllData[symbol_h][1].F,label="unaltered spectrum"); 
# scatter(Keij,rand.(Poisson.(SbgC1s+SpectrumA_1))); 
xlabel("kinetic energy [eV]",fontsize=14); ylabel("spectrum [a.u.]",fontsize=14) 
xticks(fontsize=14); yticks(fontsize=14); 
legend(fontsize=14)

ax4 = subplot(224)
title("Eph = 1216 [eV]",fontsize=14)
symbol_h = :hν_1216
plot(dictAllData[symbol_h][1].Ke,dictAllData[symbol_h][1].Sbg,label="background"); 
plot(dictAllData[symbol_h][1].Ke,dictAllData[symbol_h][1].Sbg+dictAllData[symbol_h][1].SpectrumA_1,label="noise free spectrum"); 
scatter(dictAllData[symbol_h][1].Ke,dictAllData[symbol_h][1].Snoisy,label="noisy spectrum"); 
# GeomGain = dictAllGeom[Symbol("λe_1.3000000000000003")][1].H'*dictAllGeom[Symbol("λe_1.3000000000000003")][1].ρ
# plot(dictAllData[symbol_h][1].Ke,dictAllData[symbol_h][1].Sbg+GeomGain*dictAllData[symbol_h][1].σ_cs_dens.*dictAllData[symbol_h][1].σ_tot.*dictAllData[symbol_h][1].T.*dictAllData[symbol_h][1].F,label="unaltered spectrum"); 
# scatter(Keij,rand.(Poisson.(SbgC1s+SpectrumA_1))); 
xlabel("kinetic energy [eV]",fontsize=14); ylabel("spectrum [a.u.]",fontsize=14) 
xticks(fontsize=14); yticks(fontsize=14); 
legend(fontsize=14)

tight_layout(pad=1.0, w_pad=0.5, h_pad=0.2)

# savefig("full_measurement_model.png")
# savefig("full_measurement_model.pdf")

