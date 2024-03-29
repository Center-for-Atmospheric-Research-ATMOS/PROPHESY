#------------------------------------------------------------------------------
#
# This file is part of the ATTIRE module which is licensed under CC-BY 4.0 license.
#
# Copyright (C) 2022,  Matthew Ozon, Konstantin Tumashevich and Nønne L. Prisle.
#
#------------------------------------------------------------------------------
"""
    model the spectral spread of the exciting source

    sourceSpread(hν::Cdouble,hνj::Cdouble,Δνj::Cdouble,Fνj::Cdouble)

    hνj: central photon energy
    Δνj: spectral spread of the photon beam
    Fνj: total photon flux density, e.g. m^{-2} s^{-1} if integrated over the spatial extent of the beam
"""
function sourceSpread(hν::Cdouble,hνj::Cdouble,Δνj::Cdouble,Fνj::Cdouble)
    (Fνj/sqrt(2π*Δνj^2))*exp(-(hν-hνj)^2/(2Δνj^2))
end