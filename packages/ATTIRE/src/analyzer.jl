
"""
Transmission function for the scienta R4000 kinetic enegy analyzer

T_r4000(Ke::Cdouble,E_pass::Cdouble)

Ke:     kinetic energy
E_pass: pass energy

formula from https://www.helmholtz-berlin.de/pubbin/igama_output?modus=datei&did=147
"""
function T_r4000(Ke::Cdouble,E_pass::Cdouble)
1 .-0.041*(Ke/E_pass) .+9.4e-4*(Ke/E_pass).^2 .-1.0e-5*(Ke/E_pass).^3 .+3.9e-8*(Ke/E_pass).^4
end

"""

    simulateSpectrum(Fνj::Cdouble,hνj::Cdouble,Δνj::Cdouble,Ki::Array{Cdouble,1},ΔKi::Cdouble,T::Cdouble,Kdisc::Array{Cdouble,1},Be0::Array{Cdouble,1},σ_cs_fg::Array{Cdouble,1},σ_bg_vec::Array{Cdouble,1})

    this function simulates the acquisition of a photo-electron spectrum with the parameters:

    (Fνj::Cdouble,hνj::Cdouble,Δνj::Cdouble): photon beam parameters

    (Ki::Array{Cdouble,1},ΔKi::Cdouble,T::Cdouble): analyzer parameters

    Kdisc::Array{Cdouble,1}: kinetic energy discretization 

    (Be0::Array{Cdouble,1},σ_cs_fg::Array{Cdouble,1},σ_bg_vec::Array{Cdouble,1}): cross section

    and then returns 
    Sj,Sj_bg: the expected spectrum for the kinetic energies Ki with the signal of interest (Sj) and the background (Sj_bg)
    Gm,Fl:    the discretization of the integrals over the kinetic energy and the photon energy respectively
"""
function simulateSpectrum(Fνj::Cdouble,hνj::Cdouble,Δνj::Cdouble,
Ki::Array{Cdouble,1},ΔKi::Cdouble,T::Cdouble,Kdisc::Array{Cdouble,1},
Be0::Array{Cdouble,1},σ_cs_fg::Array{Cdouble,1},σ_bg_vec::Array{Cdouble,1})

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
    Ndisc = length(Kdisc); # NOTE: the discretization of the kinetic enegy space does not have to coincide with the center of the channel, it can be whatever subdivision of that space
    NdensF = length(hνjl);
    # discretization of the integral: piecewise linear basis function
    Gm = dKe*ones(Cdouble,Ndisc); # Gm[1] = 0.5dKe; Gm[end] = 0.5dKe;
    Fl = Δhν*ones(Cdouble,NdensF); Fl[1] = 0.5Δhν; Fl[end] = 0.5Δhν;

    # define an interpolating tool whose nodes are Be0::Array{Cdouble,1},σ_cs_0::Array{Cdouble,1}
    σ_cs_interp = extrapolate(interpolate((Be0,), σ_cs_fg, Gridded(Linear())),Flat())
    σ_bg_interp = extrapolate(interpolate((Ki,), σ_bg_vec, Gridded(Linear())),Flat())

    # compute the "convolution" of the corss section by the spread of the light source and the spread of the kinetic energy analyzer
    Aij    = zeros(Cdouble,Ndisc,NdensF);      # discretization of the spread of the source and analyzer
    Sj     = Array{Cdouble,1}(undef,Nchannel); 
    Aij_bg = zeros(Cdouble,Ndisc,NdensF);      # discretization of the spread of the source and analyzer (for the background signal)
    Sj_bg  = Array{Cdouble,1}(undef,Nchannel);

    F_dens = sourceSpread.(hνjl,hνj,Δνj,Fνj)
    σ_val    = Array{Cdouble,2}(undef,NdensF,Ndisc)
    σ_bg_val = Array{Cdouble,2}(undef,NdensF,Ndisc)
    for l in 1:NdensF
        # interpolator for the current photon energy
        Be = hνjl[l] .- Kdisc;
        σ_val[l,:]    = σ_cs_interp(Be) # σ_cs_interp[Be] # no need to repeat that computation at each iteration of the i loop, once at first is enough
        σ_bg_val[l,:] = σ_bg_interp(Kdisc) # σ_bg_interp[Kdisc]    # background
    end
    σ_val[σ_val.<0.0] .= 0.0 ;
    for i in 1:Nchannel 
        φi_val = φi[i].(Kdisc)
        for m in 1:Ndisc
            for l in 1:NdensF
                Aij[m,l]    = φi_val[m]*F_dens[l]*σ_val[l,m]
                Aij_bg[m,l] = φi_val[m]*F_dens[l]*σ_bg_val[l,m]
            end
        end
        Sj[i] = Gm'*Aij*Fl;
        Sj_bg[i] = Gm'*Aij_bg*Fl;
    end

    # return cross section signal
    Sj,Sj_bg,Gm,Fl
end
