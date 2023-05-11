#------------------------------------------------------------------------------
#
# This file is part of the XPSsampling module which is licensed under CC-BY 4.0 license.
#
# Copyright (C) 2022,  Matthew Ozon, Konstantin Tumashevich and Nønne L. Prisle.
#
#------------------------------------------------------------------------------

"""
    corrCovariance(w::Array{Cdouble,1};cor_len::Cdouble=5.0)

    creates a matrix with elements: ``\\Gamma_{i,j} = w_i e^{\\frac{(i-j)^2}{0.5 cor\\_len^2}}``

    input: 

      - w: diagonal elements
      - cor_len: correlation length

    output:

      - ``\\Gamma``
"""
function corrCovariance(w::Array{Cdouble,1};cor_len::Cdouble=5.0)
    Nr = length(w);

    Γprior = zeros(Cdouble,Nr,Nr);

    for i in 1:Nr
        Γprior[i,i] = w[i]^2;
        for j in i+1:Nr
            Γprior[i,j] = Γprior[i,i]*exp(-(i-j)^2/(0.5*cor_len^2))
            Γprior[j,i] = Γprior[i,j]
        end
    end

    Γprior
end


"""
    smoothnessCovariance(w::Array{Cdouble,1};cor_len::Cdouble=5.0)

    computes three matrices 
      - symmetric matrix ``\\Gamma`` with elements ``\\Gamma_{i,j} = w_i e^{\\frac{(i-j)^2}{0.5 cor\\_len^2}}``
      - second order difference matrix ``D`` with elements ``D_{i,i} = 2, D_{i,i+1}=D_{i,i-1} = -1``
      - product matrix ``\\sqrt{\\left(D^T\\Gamma^{-1}D\\right)^{-1}}``
    
    Note: low to medium matrix size

    input: 

      - w: diagonal elements
      - cor_len: correlation length

    output:

      - ``\\sqrt{\\left(D^T\\Gamma^{-1}D\\right)^{-1}}``
      - ``\\Gamma``
      - ``D``
"""
function smoothnessCovariance(w::Array{Cdouble,1};cor_len::Cdouble=5.0)
    Nr = length(w);

    Γprior = zeros(Cdouble,Nr,Nr);
    # Dprior = D2nd(Nr+2)[:,2:end-1];
    Dprior = diagm(Nr,Nr+2,1 => 2ones(Cdouble,Nr), 0 => -ones(Cdouble,Nr) ,2 => -ones(Cdouble,Nr))[:,2:end-1]

    for i in 1:Nr
        Γprior[i,i] = w[i]^2;
        for j in i+1:Nr
            Γprior[i,j] = Γprior[i,i]*exp(-(i-j)^2/(0.5*cor_len^2))
            Γprior[j,i] = Γprior[i,j]
        end
    end

    Bprior = Dprior'*inv(Γprior)*Dprior;
    Dsqrt = real(sqrt(inv(Bprior)));
    Dsqrt, Γprior, Dprior
end
