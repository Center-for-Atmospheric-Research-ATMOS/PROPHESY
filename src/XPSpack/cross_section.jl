# load the interpolation points from file
include("C1s.jl")
# σ_C1s_interp = extrapolate(interpolate((σ_C1s[:,1],), σ_C1s[:,2], Gridded(Linear())),Flat()) # Line()
log_x_C1s = log.(σ_C1s[:,1]);
log_y_C1s = log.(σ_C1s[:,2]);
sum_xx = sum(log_x_C1s.^2);
sum_x = sum(log_x_C1s);
sum_1 = length(log_x_C1s);
sum_y = sum(log_y_C1s);
sum_xy = sum(log_y_C1s.*log_x_C1s);
μσ_C1s_bar = inv([sum_xx sum_x; sum_x sum_1])*[sum_xy; sum_y];

function σ_C1s_exp(hν_ph::Cdouble)
   exp(μσ_C1s_bar[1]*log(hν_ph) + μσ_C1s_bar[2])
end




include("O1s.jl")
# σ_O1s_interp = extrapolate(interpolate((σ_O1s[:,1],), σ_O1s[:,2], Gridded(Linear())),Flat()) # Line
log_x_O1s = log.(σ_O1s[:,1]);
log_y_O1s = log.(σ_O1s[:,2]);
sum_xx = sum(log_x_O1s.^2);
sum_x = sum(log_x_O1s);
sum_1 = length(log_x_O1s);
sum_y = sum(log_y_O1s);
sum_xy = sum(log_y_O1s.*log_x_O1s);
μσ_O1s_bar = inv([sum_xx sum_x; sum_x sum_1])*[sum_xy; sum_y];

function σ_O1s_exp(hν_ph::Cdouble)
   exp(μσ_O1s_bar[1]*log(hν_ph) + μσ_O1s_bar[2])
end


function σ_cs_orb(ħν::Cdouble,nl::String="C1s") # for the time being, only implement C1s, but you can get other orbitales from https://vuo.elettra.eu/services/elements/WebElements.html
   val = 0.0
   if (nl=="C1s")
      val = σ_C1s_exp(ħν) # σ_C1s_interp(ħν)
   else
      val = σ_O1s_exp(ħν) # σ_O1s_interp(ħν)
   end
   val
end

"""
   σ_bg(Ke::Cdouble)

   simulated total cross section of the background of an acquisition involving inelastic electron interaction
   Ke is the kinetic energy of the incoming eletron
"""
# function σ_bg(Ke::Cdouble)
#    0.05*Ke
# end
function σ_bg(Ke::Cdouble)
   10.0exp(-Ke/1000.0)
end

"""
   σ_bg_density(Ke::Array{Cdouble,1},Ke_cutoff::Cdouble,ΔKe::Cdouble)
   spectral density of the simulated cross section of the background of an acquisition involving inelastic electron interaction
   Ke is the kinetic energy of the (secondary) (re)emitted eletron
"""
function σ_bg_density(Ke::Array{Cdouble,1},Ke_cutoff::Cdouble,ΔKe::Cdouble)
   (Ke/(2.0Ke_cutoff^2))./(1.0.+exp.((Ke.-Ke_cutoff)./ΔKe))
end

function σ_bg_lin_density(Ke::Array{Cdouble,1},Ke_cutoff::Cdouble,ΔKe::Cdouble)
   (Ke[end].-Ke)./ΔKe
end