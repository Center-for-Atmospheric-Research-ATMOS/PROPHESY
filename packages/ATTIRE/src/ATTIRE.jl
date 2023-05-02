#
# ATTIRE.jl --
#
# ATTIRE.jl is a module that aims at simulating/modelling kinetic energy analyzers
# involved in XPS measurement
#
# This module depends on the package Distribution.jl (Poisson distribution)
#
#------------------------------------------------------------------------------
#
# This file is part of the ATTIRE module which is licensed under the MIT "Expat" License:
#
# Copyright (C) 2022,  Matthew Ozon, Konstantin Tumashevich and Nønne L. Prisle.
#
#------------------------------------------------------------------------------


# This file is licensed under the MIT "Expat" License:

# Copyright (c) 2022: Matthew Ozon, Konstantin Tumashevich and Nønne L. Prisle.

# Permission is hereby granted, free of charge, to any person obtaining
# a copy of this software and associated documentation files (the
# "Software"), to deal in the Software without restriction, including
# without limitation the rights to use, copy, modify, merge, publish,
# distribute, sublicense, and/or sell copies of the Software, and to
# permit persons to whom the Software is furnished to do so, subject to
# the following conditions:

# The above copyright notice and this permission notice shall be
# included in all copies or substantial portions of the Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
# EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
# MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
# IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
# CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
# TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
# SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.


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

