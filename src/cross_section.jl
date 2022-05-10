# load the interpolation points from file
include("C1s.jl")
σ_C1s_interp = extrapolate(interpolate((σ_C1s[:,1],), σ_C1s[:,2], Gridded(Linear())),Line())
function σ_cs_orb(ħν::Cdouble,nl::String="C1s") # for the time being, only implement C1s, but you can get other orbitales from https://vuo.elettra.eu/services/elements/WebElements.html
   σ_C1s_interp(ħν)
end
