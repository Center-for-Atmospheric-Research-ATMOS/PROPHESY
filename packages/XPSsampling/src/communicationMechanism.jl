# the transition mechanism generates states that samples the proposal density, which should be fairly similar to the density we hope to sample.
# In our case, the a priori distribution in the a posteriori is fairly easy to sample and so I choose to build the transmission mechanism upon it.

#TODO: optimization SA. A communication kernel that seems like a good idea is one that would start changing only the values closest to the surface and gradually changes deeper and deeper values leaving out the surface


function transmissionMechanismSmooth(x_curr::Array{Cdouble,1},Γsqrt::Array{Cdouble,2})
    x_curr + Γsqrt*randn(Cdouble,length(x_curr))
end

function transmissionMechanismBoundaries(x_curr::Array{Cdouble,1},σB::Array{Cdouble,1})
    x_prop = x_curr;
    if (rand(Cdouble)<0.5)
        x_prop[1] = x_prop[1] + σB[1]*randn(Cdouble);
    else
        x_prop[end] = x_prop[end] + σB[2]*randn(Cdouble)
    end
    x_prop
end

"""
    the transmission mechanism samples the proposal density

    The density is taken as a truncated Gaussian distribution centered on the current state x_curr
    and whose covariance matrix is Γsqrt^2. The truncation comes in to play when negative values 
    appear in the generated state
"""
function transmissionMechanism(x_curr::Array{Cdouble,1},Γsqrt::Array{Cdouble,2},σB::Array{Cdouble,1};psmooth::Cdouble=0.99)
    if (rand(Cdouble)<psmooth)
        x_prop = transmissionMechanismSmooth(x_curr,Γsqrt)
    else
        x_prop = transmissionMechanismBoundaries(x_curr,σB)
    end
    x_prop[x_prop.<0.0] .= 0.0
    x_prop
end