# load the interpolation points from file
include("C1s.jl")
σ_C1s_interp = extrapolate(interpolate((σ_C1s[:,1],), σ_C1s[:,2], Gridded(Linear())),Line())

include("O1s.jl")
σ_O1s_interp = extrapolate(interpolate((σ_O1s[:,1],), σ_O1s[:,2], Gridded(Linear())),Line())


function σ_cs_orb(ħν::Cdouble,nl::String="C1s") # for the time being, only implement C1s, but you can get other orbitales from https://vuo.elettra.eu/services/elements/WebElements.html
   val = 0.0
   if (nl=="C1s")
      val = σ_C1s_interp(ħν)
   else
      val = σ_O1s_interp(ħν)
   end
   val
end

"""
   σ_bg(Ke::Cdouble)

   simulated total cross section of the background of an acquisition involving inelastic electron interaction
   Ke is the kinetic energy of the incoming eletron
"""
function σ_bg(Ke::Cdouble)
   0.05*Ke
end

"""
   σ_bg_density(Ke::Array{Cdouble,1},Ke_cutoff::Cdouble,ΔKe::Cdouble)
   spectral density of the simulated cross section of the background of an acquisition involving inelastic electron interaction
   Ke is the kinetic energy of the (secondary) (re)emitted eletron
"""
function σ_bg_density(Ke::Array{Cdouble,1},Ke_cutoff::Cdouble,ΔKe::Cdouble)
   (Ke/(2.0Ke_cutoff^2))./(1.0.+exp.((Ke.-Ke_cutoff)./ΔKe))
end