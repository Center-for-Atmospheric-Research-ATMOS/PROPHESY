#------------------------------------------------------------------------------
#
# This file is part of the XPSpack module which is licensed under CC-BY 4.0 license.
#
# Copyright (C) 2022,  Matthew Ozon, Konstantin Tumashevich and Nønne L. Prisle.
#
#------------------------------------------------------------------------------


# Thurmer 2013
Ke_thurmer = [25.0; 30.0; 40.0; 50.0; 70.0; 100.0; 150.0; 200.0; 300.0; 400.0; 500.0; 600.0; 700.0; 800.0; 900.0]; # in eV
λ0_thurmer = [ 0.6; 0.65;  0.7;  0.8;  1.0;   1.1;   1.2;   1.5;   1.8;   2.5;   3.0;   3.8;   4.0;   5.0;   6.2]; # in nm
λe_thurmer = extrapolate(interpolate((Ke_thurmer,), λ0_thurmer, Gridded(Linear())),Line())
# 0.1sqrt.(Ke_thurmer).*(Ke_thurmer.<200.0)+0.007Ke_thurmer.*(Ke_thurmer.>=200.0)

# electron attenuation length as a function of the kinetic energy of the electron
function λe(Ke::Cdouble)
    λe_thurmer(Ke)
end

Ke_nots = [366; 623; 1031; 1600; 363; 620.0; 1031; 1583];
λe_knot = [1.33; 1.88; 2.70; 3.79; 1.32; 1.86; 2.70; 3.75];
sum_xx = sum(Ke_nots.^2);
sum_x = sum(Ke_nots);
sum_1 = length(Ke_nots);
sum_y = sum(λe_knot);
sum_xy = sum(λe_knot.*Ke_nots);

μλ_bar = inv([sum_xx sum_x; sum_x sum_1])*[sum_xy; sum_y];

function λe_exp(Ke::Cdouble)
    μλ_bar[1]*Ke + μλ_bar[2]
end
