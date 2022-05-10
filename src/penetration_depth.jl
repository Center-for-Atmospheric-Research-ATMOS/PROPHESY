# Thurmer 2013
Ke_thurmer = [25.0; 30.0; 40.0; 50.0; 70.0; 100.0; 150.0; 200.0; 300.0; 400.0; 500.0; 600.0; 700.0; 800.0; 900.0]; # in eV
λ0_thurmer = [ 0.6; 0.65;  0.7;  0.8;  1.0;   1.1;   1.2;   1.5;   1.8;   2.5;   3.0;   3.8;   4.0;   5.0;   6.2]; # in nm
λe_thurmer = extrapolate(interpolate((Ke_thurmer,), λ0_thurmer, Gridded(Linear())),Line())
# 0.1sqrt.(Ke_thurmer).*(Ke_thurmer.<200.0)+0.007Ke_thurmer.*(Ke_thurmer.>=200.0)

# electron attenuation length as a function of the kinetic energy of the electron
function λe(Ke::Cdouble)
    λe_thurmer(Ke)
end