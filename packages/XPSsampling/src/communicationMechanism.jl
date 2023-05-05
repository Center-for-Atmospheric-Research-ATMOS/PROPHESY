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

    transmissionMechanism(x_curr::Array{Cdouble,1},Γsqrt::Array{Cdouble,2},σB::Array{Cdouble,1};psmooth::Cdouble=0.99)    

    The transmission mechanism samples the proposal density which is defined by two mechanisms:
      
      - smooth: acting over all the elements of the state vector
      - boundary: acting only on the first and last element of the state vector

    The first mechanism ``q_{\\text{smooth}}`` is described by the Gaussian distribution ``\\mathcal{N}\\left(x_{\\text{curr}},\\Gamma_{\\text{sqrt}}^2\\right)``,centered on the current state ``x_{\\text{curr}}`` and whose covariance matrix is ``\\Gamma_{\\text{sqrt}}^2``.
    The second mechanism ``q_{\\text{boundary}}`` is also decribed by Gaussian distributions acting either on the first ``x^{\\text{first}}`` or the last element ``x^{\\text{last}}``.
    The proposed samples have all their entries positive.

    The overall mechanism is 
    ```math
    \\begin{equation}
        q(x_{\\text{prop}},x_{\\text{curr}}) = p_{\\text{smooth}} q_{\\text{smooth}}(\\bullet|x_{\\text{curr}}) + (1-p_{\\text{smooth}}) (\\frac{1}{2} q_{\\text{boundary}}^{\\text{first}} + \\frac{1}{2} q_{\\text{boundary}}^{\\text{last}})
    \\end{equation}
    ```
    where ``p_{\\text{smooth}}`` is the probability of choosing the mechanism 

    input: 
    
      - ``x_{\\text{curr}}`` current state
      - ``\\Gamma_{\\text{sqrt}}`` square root matrix of the covariance of the smooth mechanism ``q_{\\text{smooth}}``
      - ``\\sigma_B`` standard deviation of the boundary mechanism ``q_{\\text{boundary}}``
      - ``p_{\text{smooth}}`` probability of the smooth mechanism

    output:

      - new state
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