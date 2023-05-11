#------------------------------------------------------------------------------
#
# This file is part of the XPSpack module which is licensed under CC-BY 4.0 license.
#
# Copyright (C) 2022,  Matthew Ozon, Konstantin Tumashevich and Nønne L. Prisle.
#
#------------------------------------------------------------------------------

# the geometry structures assumes that
#  - the photon beam is along the z axis
#  - the liquid microjet is along the y axis (or the droplet axis motion, or, if no clear axis arise naturally, any direction orthogonal to the z axis)
#  - the x axis is the remaining axis that makes Oxyz an orthogonal frame of reference
#  - the reference of axis, namely O, is taken either at the center or on the surface of the sample
# maybe but not sure: the model derived from the geomStructs assumes that the incident light is orthogonal to the surface of the sample

# finger struct
"""
    fingerGeom is a mutable structure that contains the necessary parameters to describe a 1D geometry 
    where some signal may come from, hence the name finger.

    - x0,y0,z0: location of the analyzer's apperture
    - x,y:      coordinate of the electron source at the surface of the planar sample
    - z:        array of depth into the sample

    fingerGeom(x0_::Cdouble,y0_::Cdouble,z0_::Cdouble,x_::Cdouble,y_::Cdouble,z_::Array{Cdouble,1})
    creates an object of type fingerGeom
"""
mutable struct fingerGeom
    # coordinates of the analyzer's apperture (P)
    x0::Cdouble
    y0::Cdouble
    z0::Cdouble

    # coordinate of the finger at the sample's surface
    x::Cdouble
    y::Cdouble

    # probing depths along the finger (since the sampling can be non-uniform, I decide to keep an array instead of a range)
    z::Array{Cdouble,1}

    function fingerGeom()
        new(0.0,0.0,0.0, 0.0,0.0,Array{Cdouble,1}(undef,0))
    end
    function fingerGeom(x0_::Cdouble,y0_::Cdouble,z0_::Cdouble,x_::Cdouble,y_::Cdouble,z_::Array{Cdouble,1})
        new(x0_,y0_,z0_,x_,y_,z_)
    end
    function fingerGeom(ws::fingerGeom)
        new(ws.x0,ws.y0,ws.z0,ws.x,ws.y,ws.z)
    end
end

# plane struct: assumes a rectangular cuboid discretization volume
"""
    planeGeom is a mutable structure that contains the necessary parameters to describe a 3D geometry 
    where some signal may come from (some cuboid volume underneath the surface of a planar sample)

    - x0,y0,z0: location of the analyzer's apperture
    - x,y:      coordinates of the electron source at the surface of the planar sample
    - z:        array of depth into the sample

    planeGeom(x0_::Cdouble,y0_::Cdouble,z0_::Cdouble,x_::Array{Cdouble,1},y_::Array{Cdouble,1},z_::Array{Cdouble,1})
    creates an object of type planeGeom
"""
mutable struct planeGeom
    # coordinates of the analyzer's apperture (P)
    x0::Cdouble
    y0::Cdouble
    z0::Cdouble

    # coordinates of the illuminated area at the surface of the sample (the covered area may be bigger or smaller than the sample's dimension: the derived model should deal with that with the beam-light's profile)
    x::Array{Cdouble,1}
    y::Array{Cdouble,1}

    # probing depths
    z::Array{Cdouble,1}

    function planeGeom()
        new(0.0,0.0,0.0, Array{Cdouble,1}(undef,0),Array{Cdouble,1}(undef,0),Array{Cdouble,1}(undef,0))
    end
    function planeGeom(x0_::Cdouble,y0_::Cdouble,z0_::Cdouble,x_::Array{Cdouble,1},y_::Array{Cdouble,1},z_::Array{Cdouble,1})
        new(x0_,y0_,z0_,x_,y_,z_)
    end
    function planeGeom(ws::planeGeom)
        new(ws.x0,ws.y0,ws.z0,ws.x,ws.y,ws.z)
    end
end

# cylinder struct
"""
    cylinderGeom is a mutable structure that contains the necessary parameters to describe a 3D geometry 
    where some signal may come from (some cuboid volume underneath the surface of a planar sample)

    - x0,y0,z0: location of the analyzer's apperture (Cartesian coordinate)
    - μ0:       radius of the liquid microjet
    - θ,y:      coordinates of the electron source at the surface of the cylindrical sample (cylindrical coordinates)
    - r:        array of depth into the sample (radius in the sample)

    cylinderGeom(x0_::Cdouble,y0_::Cdouble,z0_::Cdouble,μ0_::Cdouble,r_::Array{Cdouble,1},θ_::Array{Cdouble,1},y_::Array{Cdouble,1})
    creates an object of type cylinderGeom
"""
mutable struct cylinderGeom
    # coordinates of the analyzer's apperture (P)
    x0::Cdouble
    y0::Cdouble
    z0::Cdouble

    # radius of the sample
    μ0::Cdouble 

    # probing depths (distance from the symmetry axis of the sample)
    r::Array{Cdouble,1}

    # coordinates of the illuminated area at the surface of the sample in polar coordinates (the covered area may be bigger or smaller than the sample's dimension: the derived model should deal with that with the beam-light's profile)
    θ::Array{Cdouble,1}
    y::Array{Cdouble,1}

    function cylinderGeom()
        new(0.0,0.0,0.0,0.0, Array{Cdouble,1}(undef,0),Array{Cdouble,1}(undef,0),Array{Cdouble,1}(undef,0))
    end
    function cylinderGeom(x0_::Cdouble,y0_::Cdouble,z0_::Cdouble,μ0_::Cdouble,r_::Array{Cdouble,1},θ_::Array{Cdouble,1},y_::Array{Cdouble,1})
        new(x0_,y0_,z0_,μ0_,r_,θ_,y_)
    end
    function cylinderGeom(ws::cylinderGeom)
        new(ws.x0,ws.y0,ws.z0,ws.μ0,ws.r,ws.θ,ws.y)
    end
end