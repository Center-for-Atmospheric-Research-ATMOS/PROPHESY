"""

   Ψ_lin_peaks_area: peak area model for a given depth discretization Zi, and XPS setup. The function will return for each element of the dictionary XPS_peak
   a photoelectric signal and peak area model (with the uncertainty model too).

   Arguments:

      Zi: an array with the depth discretization
      XPS_peak: dictionary of XPSsetup (structure containing all the elements of the setup needed to compute the model)

   Optional arguments:

      z0:    solvent surface depth (default 0.0 [nm])
      σ_z:   solvent surface width (default 0.5 [nm])
      κ_cs:  cross section relative uncertainty (default 0.05)
      κ_eal: electron attenuation length relative uncertainty (default 0.05)
"""
# function Ψ_lin_peaks_area(Zi::Array{Cdouble,1},XPS_peak::Dict{Int64,XPSsetup};σ_z::Cdouble=0.5,z0::Cdouble=0.0,κ_cs::Cdouble=0.05,κ_eal::Cdouble=0.05)
#     # measurement model
#     H_dict     = Dict{Int64,Tuple{Array{Cdouble,2},Array{Cdouble,2}}}();
#     Apeak_dict = Dict{Int64,Tuple{Array{Cdouble,1},Array{Cdouble,1}}}();
#     for (i,peak) in XPS_peak
#        # create the measurement model and uncertainty model
#        H = Ψ_lin_peaks(Zi,peak;Nz=50,σ_z=σ_z,z0=z0,κ_cs=0.0,κ_eal=0.0)
#        _,H_std = Ψ_lin_peaks_mean_and_std(Zi,peak;Nz=10,κ_cs=κ_cs,κ_eal=κ_eal,σ_z=σ_z,z0=z0);
#        # peak area (integration along the kinetic energy dimension, using Riemann sum)
#        A = dropdims(sum(H,dims=1),dims=1)*abs(peak.Ke[2]-peak.Ke[1])
#        A_std = sqrt.(dropdims(sum(H_std.^2,dims=1),dims=1))*abs(peak.Ke[2]-peak.Ke[1])
#        # set to dictionaries
#        setindex!(H_dict,(H,H_std),i);
#        setindex!(Apeak_dict,(A,A_std),i);
#     end
#     H_dict,Apeak_dict
# end
