"""

    Channel model centered on Ki [eV], with bandwidth ΔKi [eV] and gain T
    this model is over simplistic and only approximates the overall true transfer 
    function with the least amount of parameters while still carrying the idea 
    that the device blurs the signal it intends to measure
    Ki is the kinetic energy that the device intend to measure
    ΔKi is the characteristic spectrum width around Ki that will contribute to the measurement
    T is the linear gain. It depends on the collecting lens specs, the electron multiplier efficiency and the detector efficiency

    ΔKi can be set as default value to a multiple of the analyzer resolution, e.g. n*0.05, where n is 1, 2 or another integer, and 0.05 is the distance in eV between two consecutive channel centers

    φi(Ke::Cdouble,Ki::Cdouble,ΔKi::Cdouble,T::Cdouble)

    Returns the gain of the analyzer ith channel for a kinetic energy Ke. The channel
    specs are given by the center Ki, the bandwidth ΔKi and the overall gain T
"""
function φi(Ke::Cdouble,Ki::Cdouble,ΔKi::Cdouble,T::Cdouble)
    (T/sqrt(2π*ΔKi^2))*exp(-(Ke-Ki)^2/(2ΔKi^2))
end
function φi(Ke::Array{Cdouble},Ki::Cdouble,ΔKi::Cdouble,T::Cdouble)
    (T/sqrt(2π*ΔKi^2))*exp.(-(Ke.-Ki).^2/(2ΔKi^2))
end

"""
    Φi is a function that returns an array of functions, one function per channel (see φi)

    It is assumed that the energy bandwidth is the same for all channels (it is somewhat true for a given setup, at one photon energy, with a given pass energy in a given energy range, the channels can be assumed to have the same bandwidths)
    
        Φi(Ki::Array{Cdouble,1},ΔKi::Cdouble,T::Cdouble)
    Ki is an array with all the channels center
    ΔKi is the bandwidth
    T is the overall gain
"""
function Φi(Ki::Array{Cdouble,1},ΔKi::Cdouble,T::Cdouble)
    ChannelEfficiencies = Array{Function,1}(undef,length(Ki));
    for i in 1:length(Ki)
        ChannelEfficiencies[i] = (Ke::Union{Cdouble,Array{Cdouble,1}}->φi(Ke,Ki[i],ΔKi,T))
    end
    ChannelEfficiencies
end