#
# ATTIRE.jl --
#
# ATTIRE.jl is a module that aims at simulating/modelling kinetic energy analyzers
# involved in XPS measurement
#
#------------------------------------------------------------------------------
#
# This file is part of the ATTIRE module which is licensed under CC-BY 4.0 license.
#
# Copyright (C) 2022,  Matthew Ozon, Konstantin Tumashevich and Nønne L. Prisle.
#
#------------------------------------------------------------------------------


"""
This is the [`ATTIRE`](@ref), which contains
* [`ATTIRE.T_r4000`](@ref)
* [`ATTIRE.simulateSpectrum`](@ref)
* [`ATTIRE.φi`](@ref)
* [`ATTIRE.Φi`](@ref)
* [`ATTIRE.sourceSpread`](@ref)
* [`ATTIRE.countElectrons`](@ref)
"""
module ATTIRE
    # using or import what needs to be import
    using Interpolations
    import Distributions.Poisson
    export Poisson

    # here you'll find some function to simulate the analyzer/detector found in XPS experiements
    export φi, Φi, sourceSpread, countElectrons, simulateSpectrum, T_r4000;

    # the efficiency/kernel functions of the analyzer
    include("kernel.jl") 

    # here I include the source in the instrument even though it is not a part of it.
    # However, it is an important player in the spread of the measured signal.
    # it is only the spectral spread with the total amplitude, the spatial profile is left out
    # The source should be a separate package because it involves a lot of physics that is not
    # part of the analyzer (it's out of the scope of the ATTIRE package)
    include("light.jl")

    
    # model of the measurement noise (Poisson model for CCD and CEM)
    include("noise.jl")
    

    # moved to XPSmeas (in XPSpack) because it deals with sample/source alignment
    # function alignmentParameterGaussianProfile(xc::Cdouble,yc::Cdouble)
    #     # return the ratio between the integral with the photon flux and without. 
    #     # The profile is centered in (xc,yc) and is gaussian shaped
    #     # the center (xc,yc) is the offset relative to the target's center
    # end

    ##
    ## simple analyzer model
    ##
    include("analyzer.jl")

end

