#
# XPSmeas.jl --
#
# XPSmeas.jl is a module aiming at modelling XPS measurement, e.g. simulation
# of PE spectra, and inverting data, e.g. concentration profiles.
#
#
#------------------------------------------------------------------------------
#
# This file is part of the XPSmeas module which is licensed under the MIT "Expat" License:
#
# Copyright (C) 2021-2022,  Matthew Ozon.
#
#------------------------------------------------------------------------------

#TODO: split this file


mutable struct XPSdevice

   # characteristic of the experiment
   wsEXP::XPSexp                  # the relevant characteristics of a given experiment for computing the matrix model of the measurement device

   # discretization
   Nz::Int64                      # number of discretization intervals (it means that there are Nz+1 elements in the subdivision)
   Zi::Array{Cdouble,1}           # depth discretization (interval polls)
   Zi_mid::Array{Cdouble,1}       # depth discretization mid points (center of each intervals)
   ρ_tot::Array{Cdouble,1}        # total concentration profile (sum of all concentrations)
   ρ0::Cdouble                    # bulk concentration
   gz::Array{Cdouble,1}           # total relative concentration

   # the operator modelling the measurement device
   H::Array{Cdouble,2}            # Y = H*ρ_A the deterministic part of the PE signal for a given species A is the product of H and the concentration profile ρ_A

   # default ctor (it is not really meaningful, it's more for the sake of having a default constructor)
   function XPSdevice() #
     new(XPSexp(),
         0,Array{Cdouble,1}(undef,0),Array{Cdouble,1}(undef,0),Array{Cdouble,1}(undef,0),0.0,Array{Cdouble,1}(undef,0),
         Array{Cdouble,2}(undef,0,0))
   end

   # ctor: this is physically relevant (as much as the meaning of the parameters used for the model)
   function XPSdevice(ws_EXP::XPSexp,z_min::Cdouble,z_max::Cdouble,ρ_tot_exp::Array{Cdouble,1},ρ0_exp::Cdouble)
      Nz_exp = length(ρ_tot_exp)
      # compute Zi and Zi_mid
      Zi_exp = collect(range(z_min,z_max,length=Nz_exp+1));
      Zi_mid_exp = 0.5*(Zi_exp[1:end-1]+Zi_exp[2:end]);
      gz_exp = cumsum(ρ_tot_exp.*(Zi_exp[2:end]-Zi_exp[1:end-1]));
      H_exp = zeros(Cdouble,ws_EXP.Nke,Nz_exp);
      for k in 1:ws_EXP.Nke
         for i in 1:Nz_exp
            H_exp[k,i] = ws_EXP.Fν*ws_EXP.σ_β[k]*ws_EXP.αT*(exp(-gz_exp[i]/ws_EXP.λe[k])*(Zi_exp[i+1]-Zi_exp[i])) #mid-point Riemann integration
         end
      end
      new(ws_EXP,
         Nz_exp,Zi_exp,Zi_mid_exp,ρ_tot_exp,ρ0_exp,gz_exp,
         H_exp)
   end

   # cptor
   function XPSdevice(ws::XPSdevice) #
     new(ws.wsEXP,ws.Nz,copy(ws.Zi),copy(ws.Zi_mid),copy(ws.ρ_tot),ws.ρ0,copy(ws.gz), copy(ws.H),copy(ws.λe))
   end
end





# Thurmer 2013
Ke_thurmer = [25.0; 30.0; 40.0; 50.0; 70.0; 100.0; 150.0; 200.0; 300.0; 400.0; 500.0; 600.0; 700.0; 800.0; 900.0]; # in eV
λ0_thurmer = [ 0.6; 0.65;  0.7;  0.8;  1.0;   1.1;   1.2;   1.5;   1.8;   2.5;   3.0;   3.8;   4.0;   5.0;   6.2]; # in nm
λe_thurmer = extrapolate(interpolate((Ke_thurmer,), λ0_thurmer, Gridded(Linear())),Line())
# 0.1sqrt.(Ke_thurmer).*(Ke_thurmer.<200.0)+0.007Ke_thurmer.*(Ke_thurmer.>=200.0)

# electron attenuation length as a function of the kinetic energy of the electron
function λe(Ke::Cdouble)
    λe_thurmer(Ke)
end

# load the interpolation points from file
include("C1s.jl")
σ_C1s_interp = extrapolate(interpolate((σ_C1s[:,1],), σ_C1s[:,2], Gridded(Linear())),Line())
function σ_cs_orb(ħν::Cdouble,nl::String="C1s") # for the time being, only implement C1s, but you can get other orbitales from https://vuo.elettra.eu/services/elements/WebElements.html
   σ_C1s_interp(ħν)
end

# a mopdulation of the cross section spread
function Be_distribution(Be::Cdouble)
   E0 = 290.3 # 290.860 # [eV]
   E1 = 291.9
   E2 = 293.5
   γ  = E0/500.0 # that is the part that's not sure at all (maybe it can be estimated during the profile estimation)
   γ1 = E1/500.0
   γ2 = E2/500.0
   τt = [0.667503;  0.312994;  0.0195034];
   μt = [290.299;  291.868;  293.499];
   σt = [0.438926;  0.535566;  0.381694];
   x_1 = (τt[1]/σt[1])*exp.(-0.5*((Be.-μt[1])/σt[1]).^2);
   x_2 = (τt[2]/σt[2])*exp.(-0.5*((Be.-μt[2])/σt[2]).^2);
   x_3 = (τt[3]/σt[3])*exp.(-0.5*((Be.-μt[3])/σt[3]).^2);
   x_1+x_2+x_3
end

mutable struct XPSsetup

    # setup
    nl::String                     # orbitale

    Nν::Int64                      # number of probed photon energies
    ħν::Array{Cdouble,1}           # photon energy (array of length Nν)
    Fν::Array{Cdouble,1}           # photon flux w.r.t. the photon energy (array of length Nν)
    # σν::Array{Cdouble,1}           # cross section w.r.t. the photon energy for a given orbitale (array of length Nν)
    σν::Array{Cdouble,2}           # cross section w.r.t. the photon energy for a given orbitale (array of size (Nν,Nbe))

    Nke::Int64                     # number of mean kinetic energy (should be the same as Nν)
    μKe::Array{Cdouble,1}          # mean kinetic energies
    αT::Array{Cdouble,1}           # device and experiment parameter (product of α (alignment factor) and T (transmission))

    Nbe::Int64                     # number of sampled binding energies per
    Be::Array{Cdouble,2}           # binding energies (one row per photon energy, hopefully all of the same length)

    # attenuation length
    Ke::Array{Cdouble,1}           # [ħν[i].-Be[i,:] for i in 1:Nke] array of length Nν*Nbe
    λe::Array{Cdouble,1}           # 1D array for the electron attenuation length (w.r.t. K_e)


    # default ctor (it is not really meaningful, it's more for the sake of having a default constructor)
    function XPSsetup() #
        new("C1s",0,Array{Cdouble,1}(undef,0),Array{Cdouble,1}(undef,0),Array{Cdouble,1}(undef,0,0),
            0,Array{Cdouble,1}(undef,0), Array{Cdouble,1}(undef,0),
            0,Array{Cdouble,2}(undef,0,0),
            Array{Cdouble,1}(undef,0),Array{Cdouble,1}(undef,0))  # cross section and attenuation length
    end

    # ctor: this is physically relevant (as much as the meaning of the parameters used for the model)
    function XPSsetup(ħν_exp::Array{Cdouble,1},Fν_exp::Array{Cdouble,1},μKe_exp::Array{Cdouble,1},T_exp::Array{Cdouble,1},Be_exp::Array{Cdouble,2},σν_exp::Array{Cdouble,2};α_exp::Cdouble=1.0)
        N_exp = length(ħν_exp);
        if (length(Fν_exp)!=N_exp) | (length(μKe_exp)!=N_exp) | (length(T_exp)!=N_exp)
            throw("XPSsetup: not the right length")
        end
        Nbe_exp = size(Be_exp,2);
        Kes = zeros(Cdouble,N_exp*Nbe_exp);
        [Kes[(i-1)*Nbe_exp+1:i*Nbe_exp] = ħν_exp[i].-Be_exp[i,:] for i in 1:N_exp]
        new("C1s",N_exp,ħν_exp,Fν_exp,σν_exp,
            N_exp,μKe_exp,α_exp*T_exp,
            Nbe_exp,Be_exp,
            Kes,λe.(Kes))
    end

    # cptor
    function XPSsetup(ws::XPSsetup) #
        new(ws.nl,ws.Nν,ws.ħν,ws.Fν,ws.σν,
        ws.Nke,ws.μKe,ws.αT,
        ws.Nbe,ws.Be,
        ws.Ke,ws.λe)
    end
end


# integrated total concentration from the surface
function ρ_tot_int(z::Cdouble,σ_z::Cdouble=1.0,z0::Cdouble=1.0)
    softMaxA(z-z0,σ_z)
end

function e_0(z::Cdouble,z_min::Cdouble,z_max::Cdouble)
   ((z_max-z)/(z_max-z_min))*(z>=z_min)*(z<=z_max)
end

function e_k(z::Cdouble,z_min::Cdouble,z_mid::Cdouble,z_max::Cdouble)
   max(0.0,min((z-z_min)/(z_mid-z_min),(z_max-z)/(z_max-z_mid)))
end

function e_M(z::Cdouble,z_min::Cdouble,z_max::Cdouble)
   ((z-z_min)/(z_max-z_min))*(z>=z_min)*(z<=z_max)
end


## horrible way to compute an integral
function riemann(f::Function, a::Real, b::Real, n::Int; method="right")
  if method == "right"
     xs = a .+ collect(0.0:n) * (b-a)/n
     # as = [meth(f, l, r) for (l,r) in zip(xs[1:end-1], xs[2:end])]
     as = [f(r)*(r-l) for (l,r) in zip(xs[1:end-1], xs[2:end])]
  elseif method == "left"
     # meth(f,l,r) = f(l) * (r-l)
     xs = a .+ collect(0.0:n) * (b-a)/n
     as = [f(l)*(r-l) for (l,r) in zip(xs[1:end-1], xs[2:end])]
  elseif method == "trapezoid"
     # meth(f,l,r) = (1/2) * (f(l) + f(r)) * (r-l)
     xs = a .+ collect(0.0:n) * (b-a)/n
     as = [(1.0/2.0)*(f(l) + f(r))*(r-l) for (l,r) in zip(xs[1:end-1], xs[2:end])]
  elseif method == "simpsons"
     # meth(f,l,r) = (1.0/6.0) * (f(l) + 4.0*(f((l+r)/2.0)) + f(r)) * (r-l)
     xs = a .+ collect(0.0:n) * (b-a)/n
     as = [(1.0/6.0) * (f(l) + 4.0*(f((l+r)/2.0)) + f(r)) * (r-l) for (l,r) in zip(xs[1:end-1], xs[2:end])]
 else
     throw(@sprintf "quadrature %s is not implemented" method)
  end
  sum(as)
end

# cross section for a given species in a given setup #WARNING: it depends on the setup, the energy of the electron, the augular direction, and many other parameters
function σ_exp(Ke::Cdouble,μ_ke::Cdouble,σ_ke::Cdouble)
   # (1.0/(Ke*σ_ke*sqrt(2.0π)))*exp(-(log(Ke)-μ_ke)^2/(2.0σ_ke^2))
   (1.0/(σ_ke*sqrt(2.0π)))*exp(-(Ke-μ_ke)^2/(2.0σ_ke^2))
end

function σ_exp(Ke::Cdouble,peaks::Array{Tuple{Cdouble,Cdouble},1})
   val = σ_exp(Ke,peaks[1][1],peaks[1][2])
   for i in 2:length(peaks)
      val = val + σ_exp(Ke,peaks[i][1],peaks[i][2])
   end
   val
end

function σ_exp(Ke::Array{Cdouble,1},peaks::Array{Tuple{Cdouble,Cdouble},1})
   val = σ_exp.(Ke,peaks[1][1],peaks[1][2])
   for i in 2:length(peaks)
      val = val + σ_exp.(Ke,peaks[i][1],peaks[i][2])
   end
   val
end

##
## discrete approximation of the forward model
##
function Ψ(z_min::Cdouble,z_max::Cdouble,Ke::Cdouble;Nz::Int64=10,σ_z::Cdouble=2.0,z0::Cdouble=1.0,μ_ke::Cdouble=5.0,σ_ke::Cdouble=1.0)
   # riemann((x::Cdouble->exp(-x/λe(Ke))), z_min, z_max, Nz; method="simpsons")
   riemann((x::Cdouble->σ_exp(Ke,μ_ke,σ_ke)*exp(-ρ_tot_int(x,σ_z,z0)/λe(Ke))), z_min, z_max, Nz; method="simpsons")
end

function Ψ(z_min::Cdouble,z_max::Cdouble,σ_cs::Cdouble,λeal::Cdouble;Nz::Int64=10,σ_z::Cdouble=2.0,z0::Cdouble=1.0)
   # riemann((x::Cdouble->exp(-x/λe(Ke))), z_min, z_max, Nz; method="simpsons")
   riemann((x::Cdouble->σ_cs*exp(-ρ_tot_int(x,σ_z,z0)/λeal)), z_min, z_max, Nz; method="simpsons")
end

function Ψ_lin_0(z_min::Cdouble,z_max::Cdouble,Ke::Cdouble;Nz::Int64=10,σ_z::Cdouble=2.0,z0::Cdouble=1.0,μ_ke::Cdouble=5.0,σ_ke::Cdouble=1.0)
   riemann((x::Cdouble->σ_exp(Ke,μ_ke,σ_ke)*e_0(x,z_min,z_max)*exp(-ρ_tot_int(x,σ_z,z0)/λe(Ke))), z_min, z_max, Nz; method="simpsons")
end

function Ψ_lin(z_min::Cdouble,z_mid::Cdouble,z_max::Cdouble,Ke::Cdouble;Nz::Int64=10,σ_z::Cdouble=2.0,z0::Cdouble=1.0,μ_ke::Cdouble=5.0,σ_ke::Cdouble=1.0)
   riemann((x::Cdouble->σ_exp(Ke,μ_ke,σ_ke)*e_k(x,z_min,z_mid,z_max)*exp(-ρ_tot_int(x,σ_z,z0)/λe(Ke))), z_min, z_max, Nz; method="simpsons")
end

function Ψ_lin_M(z_min::Cdouble,z_max::Cdouble,Ke::Cdouble;Nz::Int64=10,σ_z::Cdouble=2.0,z0::Cdouble=1.0,μ_ke::Cdouble=5.0,σ_ke::Cdouble=1.0)
   riemann((x::Cdouble->σ_exp(Ke,μ_ke,σ_ke)*e_M(x,z_min,z_max)*exp(-ρ_tot_int(x,σ_z,z0)/λe(Ke))), z_min, z_max, Nz; method="simpsons")
end


function Ψ_lin_0(z_min::Cdouble,z_max::Cdouble,σ_cs::Cdouble,λeal::Cdouble;Nz::Int64=10,σ_z::Cdouble=2.0,z0::Cdouble=1.0)
   riemann((x::Cdouble->σ_cs*e_0(x,z_min,z_max)*exp(-ρ_tot_int(x,σ_z,z0)/λeal)), z_min, z_max, Nz; method="simpsons")
end

function Ψ_lin(z_min::Cdouble,z_mid::Cdouble,z_max::Cdouble,σ_cs::Cdouble,λeal::Cdouble;Nz::Int64=10,σ_z::Cdouble=2.0,z0::Cdouble=1.0)
   riemann((x::Cdouble->σ_cs*e_k(x,z_min,z_mid,z_max)*exp(-ρ_tot_int(x,σ_z,z0)/λeal)), z_min, z_max, Nz; method="simpsons")
end

function Ψ_lin_M(z_min::Cdouble,z_max::Cdouble,σ_cs::Cdouble,λeal::Cdouble;Nz::Int64=10,σ_z::Cdouble=2.0,z0::Cdouble=1.0)
   riemann((x::Cdouble->σ_cs*e_M(x,z_min,z_max)*exp(-ρ_tot_int(x,σ_z,z0)/λeal)), z_min, z_max, Nz; method="simpsons")
end

function Ψ(Zi::Array{Cdouble,1},Kes::Array{Cdouble,1};Nz::Int64=10,σ_z::Cdouble=2.0,z0::Cdouble=1.0,μ_ke::Cdouble=5.0,σ_ke::Cdouble=1.0)
   Nke = length(Kes);
   N = length(Zi);
   H = Array{Cdouble}(undef,Nke,N-1);
   for j in 1:Nke
      for i in 2:N
         H[j,i-1] = Ψ(Zi[i-1],Zi[i],Kes[j];Nz=Nz,σ_z=σ_z,z0=z0,μ_ke=μ_ke,σ_ke=σ_ke)
      end
   end
   H
end

function Ψ_peaks(Zi::Array{Cdouble,1},Kes::Array{Cdouble,1},peaks::Array{Tuple{Cdouble,Cdouble},1};Nz::Int64=10,σ_z::Cdouble=2.0,z0::Cdouble=1.0)
   Nke = length(Kes);
   N = length(Zi);
   H = Array{Cdouble}(undef,Nke,N-1);
   # compute the cross section with multiple peaks and the EAL
   σ_cs = σ_exp(Kes,peaks);
   λeal = λe.(Kes);
   # go for the discrete operator
   for j in 1:Nke
      for i in 2:N
         H[j,i-1] = Ψ(Zi[i-1],Zi[i],σ_cs[j],λeal[j];Nz=Nz,σ_z=σ_z,z0=z0) # Ψ(Zi[i-1],Zi[i],Kes[j];Nz=Nz,σ_z=σ_z,z0=z0)
      end
   end
   H,λeal,σ_cs
end

function Ψ_lin(Zi::Array{Cdouble,1},Kes::Array{Cdouble,1};Nz::Int64=10,σ_z::Cdouble=2.0,z0::Cdouble=1.0,μ_ke::Cdouble=5.0,σ_ke::Cdouble=1.0)
   Nke = length(Kes);
   N = length(Zi);
   H = Array{Cdouble}(undef,Nke,N);
   for j in 1:Nke
      H[j,1] = Ψ_lin_0(Zi[1],Zi[2],Kes[j];Nz=Nz,σ_z=σ_z,z0=z0,μ_ke=μ_ke,σ_ke=σ_ke)
      for i in 2:N-1
         H[j,i] = Ψ_lin(Zi[i-1],Zi[i],Zi[i+1],Kes[j];Nz=Nz,σ_z=σ_z,z0=z0,μ_ke=μ_ke,σ_ke=σ_ke)
      end
      H[j,end] = Ψ_lin_M(Zi[end-1],Zi[end],Kes[j];Nz=Nz,σ_z=σ_z,z0=z0,μ_ke=μ_ke,σ_ke=σ_ke)
   end
   H
end


function Ψ_lin_peaks(Zi::Array{Cdouble,1},Kes::Array{Cdouble,1},peaks::Array{Tuple{Cdouble,Cdouble},1};Nz::Int64=10,σ_z::Cdouble=2.0,z0::Cdouble=1.0,κ_cs::Cdouble=0.0,κ_eal::Cdouble=0.0)
   Nke = length(Kes);
   N = length(Zi);
   H = Array{Cdouble}(undef,Nke,N);
   # compute the cross section with multiple peaks and the EAL (relative error κ_cs and κ_eal)
   σ_cs = (1.0+κ_cs)*σ_exp(Kes,peaks);
   λeal = (1.0+κ_eal)*λe.(Kes);
   # go for the discrete operator
   for j in 1:Nke
      H[j,1] = Ψ_lin_0(Zi[1],Zi[2],σ_cs[j],λeal[j];Nz=Nz,σ_z=σ_z,z0=z0)
      for i in 2:N-1
         H[j,i] = Ψ_lin(Zi[i-1],Zi[i],Zi[i+1],σ_cs[j],λeal[j];Nz=Nz,σ_z=σ_z,z0=z0)
      end
      H[j,end] = Ψ_lin_M(Zi[end-1],Zi[end],σ_cs[j],λeal[j];Nz=Nz,σ_z=σ_z,z0=z0)
   end
   H,λeal,σ_cs
end

function Ψ_lin_peaks(Zi::Array{Cdouble,1},wsXPS::XPSsetup;Nz::Int64=10,σ_z::Cdouble=2.0,z0::Cdouble=1.0,κ_cs::Cdouble=0.0,κ_eal::Cdouble=0.0)
   N = length(Zi);
   H = Array{Cdouble}(undef,wsXPS.Nν*wsXPS.Nbe,N);
   # go for the discrete operator
   for j in 1:wsXPS.Nν
      for b in 1:wsXPS.Nbe #TODO: add the weight function that stands for the contribution of each binding energy in the C1s case
         H[(j-1)*wsXPS.Nbe+b,1] = Ψ_lin_0(Zi[1],Zi[2],wsXPS.αT[j]*wsXPS.Fν[j]*(1.0+κ_cs)*wsXPS.σν[j,b],(1.0+κ_eal)*wsXPS.λe[(j-1)*wsXPS.Nbe+b];Nz=Nz,σ_z=σ_z,z0=z0)
         for i in 2:N-1
            H[(j-1)*wsXPS.Nbe+b,i] = Ψ_lin(Zi[i-1],Zi[i],Zi[i+1],wsXPS.αT[j]*wsXPS.Fν[j]*(1.0+κ_cs)*wsXPS.σν[j,b],(1.0+κ_eal)*wsXPS.λe[(j-1)*wsXPS.Nbe+b];Nz=Nz,σ_z=σ_z,z0=z0)
         end
         H[(j-1)*wsXPS.Nbe+b,end] = Ψ_lin_M(Zi[end-1],Zi[end],wsXPS.αT[j]*wsXPS.Fν[j]*(1.0+κ_cs)*wsXPS.σν[j,b],(1.0+κ_eal)*wsXPS.λe[(j-1)*wsXPS.Nbe+b];Nz=Nz,σ_z=σ_z,z0=z0)
      end
   end
   H
end



function f_ij_0(z_min::Cdouble,z_max::Cdouble,λi_min::Cdouble,λi_max::Cdouble,σ_z::Cdouble,z0::Cdouble,Nz::Int64=10,Nλ::Int64=10)
   riemann((λ::Cdouble->(1.0/(λi_max-λi_min))*riemann((x::Cdouble->e_0(x,z_min,z_max)*exp(-ρ_tot_int(x,σ_z,z0)/λ)), z_min, z_max, Nz; method="simpsons")),λi_min,λi_max,Nλ;method="simpsons")
end
function f_ij(z_min::Cdouble,z_mid::Cdouble,z_max::Cdouble,λi_min::Cdouble,λi_max::Cdouble,σ_z::Cdouble,z0::Cdouble,Nz::Int64=10,Nλ::Int64=10)
   riemann((λ::Cdouble->(1.0/(λi_max-λi_min))*riemann((x::Cdouble->e_k(x,z_min,z_mid,z_max)*exp(-ρ_tot_int(x,σ_z,z0)/λ)), z_min, z_max, Nz; method="simpsons")),λi_min,λi_max,Nλ;method="simpsons")
end
function f_ij_M(z_min::Cdouble,z_max::Cdouble,λi_min::Cdouble,λi_max::Cdouble,σ_z::Cdouble,z0::Cdouble,Nz::Int64=10,Nλ::Int64=10)
   riemann((λ::Cdouble->(1.0/(λi_max-λi_min))*riemann((x::Cdouble->e_M(x,z_min,z_max)*exp(-ρ_tot_int(x,σ_z,z0)/λ)), z_min, z_max, Nz; method="simpsons")),λi_min,λi_max,Nλ;method="simpsons")
end

function f_ij_0_sq(z_min::Cdouble,z_max::Cdouble,λi_min::Cdouble,λi_max::Cdouble,σ_z::Cdouble,z0::Cdouble,Nz::Int64=10,Nλ::Int64=10)
   riemann((λ::Cdouble->(1.0/(λi_max-λi_min))*(riemann((x::Cdouble->e_0(x,z_min,z_max)*exp(-ρ_tot_int(x,σ_z,z0)/λ)), z_min, z_max, Nz; method="simpsons"))^2),λi_min,λi_max,Nλ;method="simpsons")
end
function f_ij_sq(z_min::Cdouble,z_mid::Cdouble,z_max::Cdouble,λi_min::Cdouble,λi_max::Cdouble,σ_z::Cdouble,z0::Cdouble,Nz::Int64=10,Nλ::Int64=10)
   riemann((λ::Cdouble->(1.0/(λi_max-λi_min))*(riemann((x::Cdouble->e_k(x,z_min,z_mid,z_max)*exp(-ρ_tot_int(x,σ_z,z0)/λ)), z_min, z_max, Nz; method="simpsons"))^2),λi_min,λi_max,Nλ;method="simpsons")
end
function f_ij_M_sq(z_min::Cdouble,z_max::Cdouble,λi_min::Cdouble,λi_max::Cdouble,σ_z::Cdouble,z0::Cdouble,Nz::Int64=10,Nλ::Int64=10)
   riemann((λ::Cdouble->(1.0/(λi_max-λi_min))*(riemann((x::Cdouble->e_M(x,z_min,z_max)*exp(-ρ_tot_int(x,σ_z,z0)/λ)), z_min, z_max, Nz; method="simpsons"))^2),λi_min,λi_max,Nλ;method="simpsons")
end


# uniform distribution
function Ψ_lin_peaks_mean_and_std(Zi::Array{Cdouble,1},wsXPS::XPSsetup;Nz::Int64=10,Nλ::Int64=10,σ_z::Cdouble=2.0,z0::Cdouble=1.0,κ_cs::Cdouble=0.0,κ_eal::Cdouble=0.0)
   N = length(Zi);
   H_mean = Array{Cdouble}(undef,wsXPS.Nν*wsXPS.Nbe,N);
   H_var  = Array{Cdouble}(undef,wsXPS.Nν*wsXPS.Nbe,N);
   # go for the discrete operator
   for j in 1:wsXPS.Nν
      for b in 1:wsXPS.Nbe
         λi_min = (1.0-κ_eal)*wsXPS.λe[(j-1)*wsXPS.Nbe+b];
         λi_max = (1.0+κ_eal)*wsXPS.λe[(j-1)*wsXPS.Nbe+b];
         H_mean[(j-1)*wsXPS.Nbe+b,1] = wsXPS.αT[j]*wsXPS.Fν[j]*wsXPS.σν[j,b]*f_ij_0(Zi[1],Zi[2],λi_min,λi_max,σ_z,z0,Nz,Nλ)
         H_var[(j-1)*wsXPS.Nbe+b,1]  = (wsXPS.αT[j]*wsXPS.Fν[j]*wsXPS.σν[j,b])^2*(1.0+(κ_cs^2/3.0))*f_ij_0_sq(Zi[1],Zi[2],λi_min,λi_max,σ_z,z0,Nz,Nλ)
         for i in 2:N-1
            H_mean[(j-1)*wsXPS.Nbe+b,i] = wsXPS.αT[j]*wsXPS.Fν[j]*wsXPS.σν[j,b]*f_ij(Zi[i-1],Zi[i],Zi[i+1],λi_min,λi_max,σ_z,z0,Nz,Nλ)
            H_var[(j-1)*wsXPS.Nbe+b,i]  = (wsXPS.αT[j]*wsXPS.Fν[j]*wsXPS.σν[j,b])^2*(1.0+(κ_cs^2/3.0))*f_ij_sq(Zi[i-1],Zi[i],Zi[i+1],λi_min,λi_max,σ_z,z0,Nz,Nλ)
         end
         H_mean[(j-1)*wsXPS.Nbe+b,end] = wsXPS.αT[j]*wsXPS.Fν[j]*wsXPS.σν[j,b]*f_ij_M(Zi[end-1],Zi[end],λi_min,λi_max,σ_z,z0,Nz,Nλ)
         H_var[(j-1)*wsXPS.Nbe+b,end]  = (wsXPS.αT[j]*wsXPS.Fν[j]*wsXPS.σν[j,b])^2*(1.0+(κ_cs^2/3.0))*f_ij_M_sq(Zi[end-1],Zi[end],λi_min,λi_max,σ_z,z0,Nz,Nλ)
      end
   end
   H_mean,sqrt.(H_var-H_mean.^2)
end


##
## sampling measurement operator
##

function Ψij_distribution!(KiKe::Array{Cdouble,1},zi_min::Cdouble,zi_max,λ0::Cdouble, σ0::Cdouble;κ0::Cdouble=0.05,Nz::Int64=10,σ_z::Cdouble=2.0,z0::Cdouble=1.0)
   λplus  = (1.0+κ0)*λ0;
   λminus  = (1.0-κ0)*λ0;
   σplus  = (1.0+κ0)*σ0;
   σminus = (1.0-κ0)*σ0;
   Ns = length(KiKe);
   for i in 1:Ns #MAYBE: multithread
      λi = λ0 + 2.0*(λplus-λminus)*(rand()-0.5);
      σi = σ0 + 2.0*(σplus-σminus)*(rand()-0.5);
      # KiKe[i] = riemann((x::Cdouble->σi*exp(-ρ_tot_int(x,σ_z,z0)/λi)), zi_min, zi_max, Nz; method="simpsons")
      KiKe[i] = Ψ(zi_min,zi_max,σi,λi;Nz=Nz,σ_z=σ_z,z0=z0)
   end
end

function Ψij_distribution_lin_0!(KiKe::Array{Cdouble,1},zi_min::Cdouble,zi_max,λ0::Cdouble, σ0::Cdouble;κ0::Cdouble=0.05,Nz::Int64=10,σ_z::Cdouble=2.0,z0::Cdouble=1.0)
   λplus  = (1.0+κ0)*λ0;
   λminus  = (1.0-κ0)*λ0;
   σplus  = (1.0+κ0)*σ0;
   σminus = (1.0-κ0)*σ0;
   Ns = length(KiKe);
   for i in 1:Ns #MAYBE: multithread
      λi = λ0 + 2.0*(λplus-λminus)*(rand()-0.5);
      σi = σ0 + 2.0*(σplus-σminus)*(rand()-0.5);
      # KiKe[i] = riemann((x::Cdouble->σi*exp(-ρ_tot_int(x,σ_z,z0)/λi)), zi_min, zi_max, Nz; method="simpsons")
      KiKe[i] = Ψ_lin_0(zi_min,zi_max,σi,λi;Nz=Nz,σ_z=σ_z,z0=z0)
   end
end
function Ψij_distribution_lin!(KiKe::Array{Cdouble,1},zi_min::Cdouble,zi_mid::Cdouble,zi_max,λ0::Cdouble, σ0::Cdouble;κ0::Cdouble=0.05,Nz::Int64=10,σ_z::Cdouble=2.0,z0::Cdouble=1.0)
   λplus  = (1.0+κ0)*λ0;
   λminus  = (1.0-κ0)*λ0;
   σplus  = (1.0+κ0)*σ0;
   σminus = (1.0-κ0)*σ0;
   Ns = length(KiKe);
   for i in 1:Ns #MAYBE: multithread
      λi = λ0 + 2.0*(λplus-λminus)*(rand()-0.5);
      σi = σ0 + 2.0*(σplus-σminus)*(rand()-0.5);
      # KiKe[i] = riemann((x::Cdouble->σi*exp(-ρ_tot_int(x,σ_z,z0)/λi)), zi_min, zi_max, Nz; method="simpsons")
      KiKe[i] = Ψ_lin(zi_min,zi_mid,zi_max,σi,λi;Nz=Nz,σ_z=σ_z,z0=z0)
   end
end
function Ψij_distribution_lin_M!(KiKe::Array{Cdouble,1},zi_min::Cdouble,zi_max,λ0::Cdouble, σ0::Cdouble;κ0::Cdouble=0.05,Nz::Int64=10,σ_z::Cdouble=2.0,z0::Cdouble=1.0)
   λplus  = (1.0+κ0)*λ0;
   λminus  = (1.0-κ0)*λ0;
   σplus  = (1.0+κ0)*σ0;
   σminus = (1.0-κ0)*σ0;
   Ns = length(KiKe);
   for i in 1:Ns #MAYBE: multithread
      λi = λ0 + 2.0*(λplus-λminus)*(rand()-0.5);
      σi = σ0 + 2.0*(σplus-σminus)*(rand()-0.5);
      # KiKe[i] = riemann((x::Cdouble->σi*exp(-ρ_tot_int(x,σ_z,z0)/λi)), zi_min, zi_max, Nz; method="simpsons")
      KiKe[i] = Ψ_lin_M(zi_min,zi_max,σi,λi;Nz=Nz,σ_z=σ_z,z0=z0)
   end
end

# Ns: number of samples for the parameters for each matrix elements
# Zi: discretization of the depth (the boundaries of the intervals used for the integrations because of the rectangle functions)
# Kes: kinetic energy discretization (mid-points)
# Nz: number of points used in the Simpsons quadrature method
# z0,σ_z: total concentration profil characteristics
# μ_ke,σ_ke: cross section characteristics
function mean_and_std_map(Ns::Int64,Zi::Array{Cdouble,1},Kes::Array{Cdouble,1};Nz::Int64=10,σ_z::Cdouble=2.0,z0::Cdouble=1.0,μ_ke::Cdouble=5.0,σ_ke::Cdouble=2.0,κ0::Cdouble=0.05)
   Nke = length(Kes);
   Ni  = length(Zi);
   H_mean = zeros(Cdouble,Nke,Ni-1);
   H_std  = zeros(Cdouble,Nke,Ni-1);
   KiKe = zeros(Cdouble,Ns);
   for i in 1:Nke  #MAYBE: multithread
      λ0     = λe(Kes[i]);
      σ0     = σ_exp(Kes[i],μ_ke,σ_ke);
      for j in 1:Ni-1
         # draw samples
         Ψij_distribution!(KiKe,Zi[j], Zi[j+1],λ0, σ0;κ0=κ0,Nz=Nz,σ_z=σ_z,z0=z0)
         # compute mean and std
         H_mean[i,j] = mean(KiKe);
         H_std[i,j]  = std(KiKe);
      end
   end
   # return mean and std
   H_mean,H_std
end

function mean_and_std_map_peaks(Ns::Int64,Zi::Array{Cdouble,1},Kes::Array{Cdouble,1},peaks::Array{Tuple{Cdouble,Cdouble},1};Nz::Int64=10,σ_z::Cdouble=2.0,z0::Cdouble=1.0,κ0::Cdouble=0.05)
   Nke = length(Kes);
   Ni  = length(Zi);
   H_mean = zeros(Cdouble,Nke,Ni-1);
   H_std  = zeros(Cdouble,Nke,Ni-1);
   KiKe = zeros(Cdouble,Ns);
   # compute the cross section with multiple peaks and the EAL
   σ_cs = σ_exp(Kes,peaks);
   λeal = λe.(Kes);
   for i in 1:Nke  #MAYBE: multithread
      λ0     = λeal[i]; # λe(Kes[i]);
      σ0     = σ_cs[i]; # σ_exp(Kes[i],μ_ke,σ_ke);
      for j in 1:Ni-1
         # draw samples
         Ψij_distribution!(KiKe,Zi[j], Zi[j+1],λ0, σ0;κ0=κ0,Nz=Nz,σ_z=σ_z,z0=z0)
         # compute mean and std
         H_mean[i,j] = mean(KiKe);
         H_std[i,j]  = std(KiKe);
      end
   end
   # return mean and std
   H_mean,H_std
end

function mean_and_std_map_lin(Ns::Int64,Zi::Array{Cdouble,1},Kes::Array{Cdouble,1};Nz::Int64=10,σ_z::Cdouble=2.0,z0::Cdouble=1.0,μ_ke::Cdouble=5.0,σ_ke::Cdouble=2.0,κ0::Cdouble=0.05)
   Nke = length(Kes);
   Ni  = length(Zi);
   H_mean = zeros(Cdouble,Nke,Ni);
   H_std  = zeros(Cdouble,Nke,Ni);
   KiKe = zeros(Cdouble,Ns);
   for i in 1:Nke  #MAYBE: multithread
      λ0     = λe(Kes[i]);
      σ0     = σ_exp(Kes[i],μ_ke,σ_ke);
      # draw samples
      Ψij_distribution_lin_0!(KiKe,Zi[1],Zi[2],λ0, σ0;κ0=κ0,Nz=Nz,σ_z=σ_z,z0=z0)
      # compute mean and std
      H_mean[i,1] = mean(KiKe);
      H_std[i,1]  = std(KiKe);
      for j in 2:Ni-1
         # draw samples
         Ψij_distribution_lin!(KiKe,Zi[j-1],Zi[j],Zi[j+1],λ0, σ0;κ0=κ0,Nz=Nz,σ_z=σ_z,z0=z0)
         # compute mean and std
         H_mean[i,j] = mean(KiKe);
         H_std[i,j]  = std(KiKe);
      end
      # draw samples
      Ψij_distribution_lin_M!(KiKe,Zi[end-1],Zi[end],λ0, σ0;κ0=κ0,Nz=Nz,σ_z=σ_z,z0=z0)
      # compute mean and std
      H_mean[i,end] = mean(KiKe);
      H_std[i,end]  = std(KiKe);
   end
   # return mean and std
   H_mean,H_std
end

function mean_and_std_map_lin_peaks(Ns::Int64,Zi::Array{Cdouble,1},Kes::Array{Cdouble,1},peaks::Array{Tuple{Cdouble,Cdouble},1};Nz::Int64=10,σ_z::Cdouble=2.0,z0::Cdouble=1.0,κ0::Cdouble=0.05)
   Nke = length(Kes);
   Ni  = length(Zi);
   H_mean = zeros(Cdouble,Nke,Ni);
   H_std  = zeros(Cdouble,Nke,Ni);
   KiKe = zeros(Cdouble,Ns);
   # compute the cross section with multiple peaks and the EAL
   σ_cs = σ_exp(Kes,peaks);
   λeal = λe.(Kes);
   for i in 1:Nke  #MAYBE: multithread
      λ0     = λeal[i]; # λe(Kes[i]);
      σ0     = σ_cs[i]; # σ_exp(Kes[i],μ_ke,σ_ke);
      # draw samples
      Ψij_distribution_lin_0!(KiKe,Zi[1],Zi[2],λ0, σ0;κ0=κ0,Nz=Nz,σ_z=σ_z,z0=z0)
      # compute mean and std
      H_mean[i,1] = mean(KiKe);
      H_std[i,1]  = std(KiKe);
      for j in 2:Ni-1
         # draw samples
         Ψij_distribution_lin!(KiKe,Zi[j-1],Zi[j],Zi[j+1],λ0, σ0;κ0=κ0,Nz=Nz,σ_z=σ_z,z0=z0)
         # compute mean and std
         H_mean[i,j] = mean(KiKe);
         H_std[i,j]  = std(KiKe);
      end
      # draw samples
      Ψij_distribution_lin_M!(KiKe,Zi[end-1],Zi[end],λ0, σ0;κ0=κ0,Nz=Nz,σ_z=σ_z,z0=z0)
      # compute mean and std
      H_mean[i,end] = mean(KiKe);
      H_std[i,end]  = std(KiKe);
   end
   # return mean and std
   H_mean,H_std
end

# error due to errors in the EAL #WARNING: deprecated
function ΔΨλ(Zi::Array{Cdouble,1},Kes::Array{Cdouble,1},ΔEALs::Array{Cdouble,1};Nz::Int64=10,σ_z::Cdouble=2.0,z0::Cdouble=1.0,μ_ke::Cdouble=5.0,σ_ke::Cdouble=2.0)
   Nke = length(Kes);
   N = length(Zi);
   H = Array{Cdouble}(undef,Nke,N-1);
   for j in 1:Nke
      for i in 2:N
         H[j,i-1] = riemann((x::Cdouble->(ΔEALs[j]/λe(Kes[j]))*(ρ_tot_int(x,σ_z,z0)/λe(Kes[j]))*σ_exp(Kes[j],μ_ke,σ_ke)*exp(-ρ_tot_int(x,σ_z,z0)/λe(Kes[j]))), Zi[i-1], Zi[i], Nz; method="simpsons")
      end
   end
   H
end

# error due to errors in the cross section #WARNING: deprecated
function ΔΨσ(Zi::Array{Cdouble,1},Kes::Array{Cdouble,1},Δσ::Array{Cdouble,1};Nz::Int64=10,σ_z::Cdouble=2.0,z0::Cdouble=1.0,μ_ke::Cdouble=5.0,σ_ke::Cdouble=2.0)
   Nke = length(Kes);
   N = length(Zi);
   H = Array{Cdouble}(undef,Nke,N-1);
   for j in 1:Nke
      for i in 2:N
         H[j,i-1] = riemann((x::Cdouble->Δσ[j]*exp(-ρ_tot_int(x,σ_z,z0)/λe(Kes[j]))), Zi[i-1],Zi[i], Nz; method="simpsons")
      end
   end
   H
end


##
## discretization uncertainty #WARNING: deprecated
##


function δintegrals(Zi::Array{Cdouble,1},Kes::Array{Cdouble,1};Nz::Int64=10,σ_z::Cdouble=2.0,z0::Cdouble=1.0,μ_ke::Cdouble=5.0,σ_ke::Cdouble=2.0)
   Nke = length(Kes);
   N = length(Zi);
   integral_0 = zeros(Cdouble,Nke,N-1);
   integral_1 = zeros(Cdouble,Nke,N-1);
   for i in 1:Nke
      for j in 1:(N-1)
         z_mid = 0.5*(Zi[j]+Zi[j+1]);
         integral_0[i,j] = riemann((x::Cdouble->(x-z_mid)*σ_exp(Kes[i],μ_ke,σ_ke)*exp(-ρ_tot_int(x,σ_z,z0)/λe(Kes[i]))), Zi[j], Zi[j+1], Nz; method="simpsons");
         integral_1[i,j] = riemann((x::Cdouble->0.5*((x-z_mid)^2)*σ_exp(Kes[i],μ_ke,σ_ke)*exp(-ρ_tot_int(x,σ_z,z0)/λe(Kes[i]))), Zi[j], Zi[j+1], Nz; method="simpsons");
      end
   end
   integral_0,integral_1
end

function δintegrals(Zi::Array{Cdouble,1},Kes::Array{Cdouble,1},peaks::Array{Tuple{Cdouble,Cdouble},1};Nz::Int64=10,σ_z::Cdouble=2.0,z0::Cdouble=1.0)
   Nke = length(Kes);
   N = length(Zi);
   integral_0 = zeros(Cdouble,Nke,N-1);
   integral_1 = zeros(Cdouble,Nke,N-1);
   σ_cs = σ_exp(Kes,peaks);
   λeal = λe.(Kes)
   for i in 1:Nke
      for j in 1:(N-1)
         z_mid = 0.5*(Zi[j]+Zi[j+1]);
         integral_0[i,j] = riemann((x::Cdouble->(x-z_mid)*σ_cs[i]*exp(-ρ_tot_int(x,σ_z,z0)/λeal[i])), Zi[j], Zi[j+1], Nz; method="simpsons");
         integral_1[i,j] = riemann((x::Cdouble->0.5*((x-z_mid)^2)*σ_cs[i]*exp(-ρ_tot_int(x,σ_z,z0)/λeal[i])), Zi[j], Zi[j+1], Nz; method="simpsons");
      end
   end
   integral_0,integral_1
end
