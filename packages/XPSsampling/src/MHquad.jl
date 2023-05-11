#------------------------------------------------------------------------------
#
# This file is part of the XPSsampling module which is licensed under CC-BY 4.0 license.
#
# Copyright (C) 2022,  Matthew Ozon, Konstantin Tumashevich and Nønne L. Prisle.
#
#------------------------------------------------------------------------------

"""
    acceptSample(ρ_cur::Array{Cdouble,1},ρ_prop::Array{Cdouble,1},y::Array{Cdouble,1},yd::Array{Cdouble,1},ΓIinv::Array{Cdouble,2},Γdinv::Array{Cdouble,2},H::Array{Cdouble,2},D::Array{Cdouble,2})

    Checks out the proposed state ``\\rho_{\\text{prop}}`` compared with the current state ``\\rho_{\\text{curr}}`` in terms of the energy function ``U`` to decide which state to accept.
    The proposed state is accepted with probability ``\\masthsf{P}(\\rho_{\\text{prop}}) = \\min(1,e^{-\\Delta})`` where ``\\Delta = U(\\rho_{\\text{prop}}) - U(\\rho_{\\text{curr}})``

    The cost function ``U`` is 
    ```math
        \\begin{equation}
            U(\\rho) = ||y-H\\rho||_{\\Gamma_I}^2 + ||y_d-D\\rho||_{\\Gamma_D}^2
        \\end{equation}
    ```
    where ``\\|y-H\\rho\\|_{\\Gamma_I}^2`` is the data fidelity term with ``y`` the vector of measurement data, ``H`` the measurement operator and ``\\Gamma_I`` the noise covariance matrix.
    The term ``\\|y_d-D\\rho\\|_{\\Gamma_D}^2`` is the regularizer of the inverse problem associated to the cost function ``U`` with ``D`` the regularization operator (e.g. second order difference), ``y_d`` a vector of expected values of the regularization and ``\\Gamma_D`` the covariance matrix of the regularization (strength and correlation)

    input:

      - ``\\rho_{\\test{cur}}`` and ``\\rho_{\\test{prop}}``: current and proposed states
      - ``y``: measurement vector
      - ``y_d``: regularization expected values
      - ``\\Gamma_I^{-1}``: inverse covariance matrix of the measurement
      - ``\\Gamma_D^{-1}``: inverse covariance matrix of the regularization
      - ``H``: measurement operator
      - ``D``: regularization operator

    output: 

      - ``\\rho_{\\text{new}}``: the new state
      - ``-\\Delta``: the cost function variation
"""
function acceptSample(ρ_cur::Array{Cdouble,1},ρ_prop::Array{Cdouble,1},y::Array{Cdouble,1},yd::Array{Cdouble,1},ΓIinv::Array{Cdouble,2},Γdinv::Array{Cdouble,2},H::Array{Cdouble,2},D::Array{Cdouble,2})
    # if the posterior probability is larger for the proposed state ρ_prop than the current state ρ_cur, then accept the state, otherwise, reject it with probability p

    # likelihood
    r_cp = 0.5*(y-H*ρ_cur)'*ΓIinv*(y-H*ρ_cur)-0.5*(y-H*ρ_prop)'*ΓIinv*(y-H*ρ_prop);
    # smoothness
    r_cp = r_cp + 0.5*(yd-D*ρ_cur)'*Γdinv*(yd-D*ρ_cur)-0.5*(yd-D*ρ_prop)'*Γdinv*(yd-D*ρ_prop); 
    
    ρ_new = ρ_cur;
    if (r_cp>=0.0)
        # unconditionally accept the new state
        ρ_new = ρ_prop
    else
        # accept the state with probability e^{r_cp}
        if (log(rand())<=r_cp)
            ρ_new = ρ_prop
        else
            r_cp = 0.0
        end
    end
    ρ_new,r_cp 
end


"""
    samplePosterior(ρ_start::Array{Cdouble,1},Γsqrt::Array{Cdouble,2},y::Array{Cdouble,1},yd::Array{Cdouble,1},ΓIinv::Array{Cdouble,2},Γdinv::Array{Cdouble,2},H::Array{Cdouble,2},D::Array{Cdouble,2};Ns::Int64=10000)

    computes a series of samples ``(\\rho_i)_i`` drawn from the generative communication mechanism ``q(\\rho_{\\text{prop}},\\rho_{\\text{curr}})`` and the cost function variations ``(\\Delta_{i})_{i}`` associated with each move.
    The cost function is 
    ```math
        \\begin{equation}
            U(\\rho) = ||y-H\\rho||_{\\Gamma_I}^2 + ||y_d-D\\rho||_{\\Gamma_D}^2
        \\end{equation}
    ```
    where ``\\|y-H\rho\\|_{\\Gamma_I}^2`` is the data fidelity term with ``y`` the vector of measurement data, ``H`` the measurement operator and ``\\Gamma_I`` the noise covariance matrix.
    The term ``\\|y_d-D\rho\\|_{\\Gamma_D}^2`` is the regularizer of the inverse problem associated to the cost function ``U`` with ``D`` the regularization operator (e.g. second order difference), ``y_d`` a vector of expected values of the regularization and ``\\Gamma_D`` the covariance matrix of the regularization (strength and correlation)

    At each step, a state is proposed ``\\rho_{\\text{prop}} \\sim q(\\rho_{\\bullet,\\rho_{\\text{curr}})`` and the state is accepted with probability ``\\masthsf{P}(\\rho_i) = \\min(1,e^{-\\Delta_{i}})``

    input:

      - ``\\rho_{\\test{start}}``: initial state of the chain of states
      - ``\\Gamma_{\\text{sqrt}}``: square root matrix of the covariance of the smooth mechanism ``q_{\\text{smooth}}``
      - ``y``: measurement vector
      - ``y_d``: regularization expected values
      - ``\\Gamma_I^{-1}``: inverse covariance matrix of the measurement
      - ``\\Gamma_D^{-1}``: inverse covariance matrix of the regularization
      - ``H``: measurement operator
      - ``D``: regularization operator
      - ``N_s``: chain length

    output: 

      - ``(\\rho_i)_i``: chain of states
      - ``(\\Delta_{i})_{i}``: cost function variation associated to the chain of states
"""
function samplePosterior(ρ_start::Array{Cdouble,1},Γsqrt::Array{Cdouble,2},y::Array{Cdouble,1},yd::Array{Cdouble,1},ΓIinv::Array{Cdouble,2},Γdinv::Array{Cdouble,2},H::Array{Cdouble,2},D::Array{Cdouble,2};Ns::Int64=10000)
    # all samples
    ρ_all = zeros(Cdouble,Ns+1,length(ρ_start))
    ρ_all[1,:] = ρ_start;
    deltaU = zeros(Cdouble,Ns);
    for i in 1:Ns
        # draw a new sample from a distribution not to far from the actual one
        ρ_all[i+1,:] = transmissionMechanismSmooth(ρ_all[i,:],Γsqrt)
        ρ_all[i+1,ρ_all[i+1,:].<0.0] .= 0.0  
        
        # accept or reject the sample
        ρ_all[i+1,:],deltaU[i] = acceptSample(ρ_all[i,:],ρ_all[i+1,:],y,yd,ΓIinv,Γdinv,H,D)
    end
    ρ_all,deltaU
end

"""
    samplePosteriorMeanAndCov(ρ_start::Array{Cdouble,1},Γsqrt::Array{Cdouble,2},y::Array{Cdouble,1},yd::Array{Cdouble,1},ΓIinv::Array{Cdouble,2},Γdinv::Array{Cdouble,2},H::Array{Cdouble,2},D::Array{Cdouble,2};Ns::Int64=10000,Nburn::Int64=1000)

    computes the mean and covariance matrix from a series of samples ``(\\rho_i)_i`` drawn from the generative communication mechanism ``q(\\rho_{\\text{prop}},\\rho_{\\text{curr}})`` and the cost function variations ``(\\Delta_{i})_{i}`` associated with each move.
    The cost function is 
    ```math
        \\begin{equation}
            U(\\rho) = ||y-H\\rho||_{\\Gamma_I}^2 + ||y_d-D\\rho||_{\\Gamma_D}^2
        \\end{equation}
    ```
    where ``\\|y-H\rho\\|_{\\Gamma_I}^2`` is the data fidelity term with ``y`` the vector of measurement data, ``H`` the measurement operator and ``\\Gamma_I`` the noise covariance matrix.
    The term ``\\|y_d-D\rho\\|_{\\Gamma_D}^2`` is the regularizer of the inverse problem associated to the cost function ``U`` with ``D`` the regularization operator (e.g. second order difference), ``y_d`` a vector of expected values of the regularization and ``\\Gamma_D`` the covariance matrix of the regularization (strength and correlation)

    At each step, a state is proposed ``\\rho_{\\text{prop}} \\sim q(\\rho_{\\bullet,\\rho_{\\text{curr}})`` and the state is accepted with probability ``\\masthsf{P}(\\rho_i) = \\min(1,e^{-\\Delta_{i}})``

    input:

      - ``\\rho_{\\test{start}}``: initial state of the chain of states
      - ``\\Gamma_{\\text{sqrt}}``: square root matrix of the covariance of the smooth mechanism ``q_{\\text{smooth}}``
      - ``y``: measurement vector
      - ``y_d``: regularization expected values
      - ``\\Gamma_I^{-1}``: inverse covariance matrix of the measurement
      - ``\\Gamma_D^{-1}``: inverse covariance matrix of the regularization
      - ``H``: measurement operator
      - ``D``: regularization operator
      - ``N_s``: chain length
      - ``N_{\\text{burn}}``: burn in period length

    output: 

      - ``\\mu_{\\rho|H,y}`` and ``\\Gamma_{\\rho|H,y}``: estimated mean state and covariance matrix from the samples `(\\rho_i)_{N_{\\text{burn}}<i<N_{s}}``
      - ``(\\Delta_{i})_{i}``: cost function variation associated to the chain of states ``(\\rho_i)_i``
"""
function samplePosteriorMeanAndCov(ρ_start::Array{Cdouble,1},Γsqrt::Array{Cdouble,2},y::Array{Cdouble,1},yd::Array{Cdouble,1},ΓIinv::Array{Cdouble,2},Γdinv::Array{Cdouble,2},H::Array{Cdouble,2},D::Array{Cdouble,2};Ns::Int64=10000,Nburn::Int64=1000)
    # all samples
    ρ_curr = copy(ρ_start);
    ρ_prop = copy(ρ_start);
    μρ_cum = zeros(Cdouble,length(ρ_start));
    Γρ_cum = zeros(Cdouble,length(ρ_start),length(ρ_start));
    deltaU = zeros(Cdouble,Ns);
    for i in 1:Ns
        # global ρ_prop
        # global deltaU
        # draw a new sample from a distribution not to far from the actual one
        ρ_prop = transmissionMechanismSmooth(ρ_curr,Γsqrt)
        ρ_prop[ρ_prop.<0.0] .= 0.0  
        
        # accept or reject the sample
        ρ_prop,deltaU[i] = acceptSample(ρ_curr,ρ_prop,y,yd,ΓIinv,Γdinv,H,D)
        if (i>=Nburn)
            μρ_cum = μρ_cum + ρ_prop; # global 
            Γρ_cum = Γρ_cum + ρ_prop*ρ_prop'; #global 
        end
        ρ_curr = ρ_prop;
    end
    μρ_cum = μρ_cum/(Ns-Nburn+1);
    Γρ_cum = (1.0/(Ns-Nburn))*(Γρ_cum - (Ns-Nburn+1)*μρ_cum*μρ_cum');
    μρ_cum,Γρ_cum,deltaU
end



"""
    acceptSampleMargin(ρ_cur::Array{Cdouble,1},ρ_prop::Array{Cdouble,1},y::Array{Cdouble,1},yd::Array{Cdouble,1},ΓIinv::Array{Cdouble,2},Γdinv::Array{Cdouble,2},H::Array{Cdouble,2},D::Array{Cdouble,2},ΓHΓyinv::Array{Cdouble,2})

    Checks out the proposed state ``\\rho_{\\text{prop}}`` compared with the current state ``\\rho_{\\text{curr}}`` in terms of the energy function ``U`` to decide which state to accept.
    The proposed state is accepted with probability ``\\masthsf{P}(\\rho_{\\text{prop}}) = \\min(1,e^{-\\Delta})`` where ``\\Delta = U(\\rho_{\\text{prop}}) - U(\\rho_{\\text{curr}})``

    The cost function ``U`` is 
    ```math
        \\begin{equation}
            U(\\rho) = ||y-H\\rho||_{\\Gamma_I}^2 + ||y_d-D\\rho||_{\\Gamma_D}^2 + \\rho^T\\Gamma_{H,y}\\rho
        \\end{equation}
    ```
    where ``\\|y-H\\rho\\|_{\\Gamma_I}^2`` is the data fidelity term with ``y`` the vector of measurement data, ``H`` the measurement operator and ``\\Gamma_I`` the noise covariance matrix.
    The term ``\\|y_d-D\\rho\\|_{\\Gamma_D}^2`` is the regularizer of the inverse problem associated to the cost function ``U`` with ``D`` the regularization operator (e.g. second order difference), ``y_d`` a vector of expected values of the regularization and ``\\Gamma_D`` the covariance matrix of the regularization (strength and correlation)
    The last term is related to the modelling uncertainties, i.e. errors in the measurement operator. The covaraince matrix ``\\Gamma_{H,y}`` is defined as the sum of the contribution of the error for each measurement normalized by the noise variance.
    Formaly, 
    ```math
        \\begin{equation}
            \\Gamma_{H,y} = \\underset{k}{\\sum} \\frac{1}{\\sigma_k^2} \\Gamma_{H_k}
        \\end{equation}
    ```
    where ``\\sigma_k^2`` is the variance of the ``k^{\\text{th}}`` measurement and ``\\Gamma_{H_k}`` is the covariance matrix of the measurement model the same measurement.

    input:

      - ``\\rho_{\\test{cur}}`` and ``\\rho_{\\test{prop}}``: current and proposed states
      - ``y``: measurement vector
      - ``y_d``: regularization expected values
      - ``\\Gamma_I^{-1}``: inverse covariance matrix of the measurement
      - ``\\Gamma_D^{-1}``: inverse covariance matrix of the regularization
      - ``H``: measurement operator
      - ``D``: regularization operator
      - ``\\Gamma_{H,y}``: operator-to-noise covariance matrix

    output: 

      - ``\\rho_{\\text{new}}``: the new state
      - ``-\\Delta``: the cost function variation
"""
function acceptSampleMargin(ρ_cur::Array{Cdouble,1},ρ_prop::Array{Cdouble,1},y::Array{Cdouble,1},yd::Array{Cdouble,1},ΓIinv::Array{Cdouble,2},Γdinv::Array{Cdouble,2},H::Array{Cdouble,2},D::Array{Cdouble,2},ΓHΓyinv::Array{Cdouble,2})
    # if the posterior probability is larger for the proposed state ρ_prop than the current state ρ_cur, then accept the state, otherwise, reject it with probability p

    # likelihood
    r_cp = 0.5*(y-H*ρ_cur)'*ΓIinv*(y-H*ρ_cur)-0.5*(y-H*ρ_prop)'*ΓIinv*(y-H*ρ_prop);
    # smoothness
    r_cp = r_cp + 0.5*(yd-D*ρ_cur)'*Γdinv*(yd-D*ρ_cur)-0.5*(yd-D*ρ_prop)'*Γdinv*(yd-D*ρ_prop); 
    # marginalization of errors
    r_cp = r_cp + 0.5*ρ_cur'*ΓHΓyinv*ρ_cur-0.5*ρ_prop'*ΓHΓyinv*ρ_prop; 
    
    ρ_new = ρ_cur;
    if (r_cp>=0.0)
        # unconditionally accept the new state
        ρ_new = ρ_prop
    else
        # accept the state with probability e^{r_cp}
        if (log(rand())<=r_cp)
            ρ_new = ρ_prop
        else
            r_cp = 0.0
        end
    end
    ρ_new,r_cp # maybe it could return the computed values
end



"""
    samplePosteriorMargin(ρ_start::Array{Cdouble,1},Γsqrt::Array{Cdouble,2},y::Array{Cdouble,1},yd::Array{Cdouble,1},ΓIinv::Array{Cdouble,2},Γdinv::Array{Cdouble,2},H::Array{Cdouble,2},D::Array{Cdouble,2},ΓHΓyinv::Array{Cdouble,2};Ns::Int64=10000)

    computes a series of samples ``(\\rho_i)_i`` drawn from the generative communication mechanism ``q(\\rho_{\\text{prop}},\\rho_{\\text{curr}})`` and the cost function variations ``(\\Delta_{i})_{i}`` associated with each move.
    The cost function is 
    ```math
        \\begin{equation}
            U(\\rho) = ||y-H\\rho||_{\\Gamma_I}^2 + ||y_d-D\\rho||_{\\Gamma_D}^2
        \\end{equation}
    ```
    where ``\\|y-H\\rho\\|_{\\Gamma_I}^2`` is the data fidelity term with ``y`` the vector of measurement data, ``H`` the measurement operator and ``\\Gamma_I`` the noise covariance matrix.
    The term ``\\|y_d-D\\rho\\|_{\\Gamma_D}^2`` is the regularizer of the inverse problem associated to the cost function ``U`` with ``D`` the regularization operator (e.g. second order difference), ``y_d`` a vector of expected values of the regularization and ``\\Gamma_D`` the covariance matrix of the regularization (strength and correlation)
    The last term is related to the modelling uncertainties, i.e. errors in the measurement operator. The covaraince matrix ``\\Gamma_{H,y}`` is defined as the sum of the contribution of the error for each measurement normalized by the noise variance.
    Formaly, 
    ```math
        \\begin{equation}
            \\Gamma_{H,y} = \\underset{k}{\\sum} \\frac{1}{\\sigma_k^2} \\Gamma_{H_k}
        \\end{equation}
    ```
    where ``\\sigma_k^2`` is the variance of the ``k^{\\text{th}}`` measurement and ``\\Gamma_{H_k}`` is the covariance matrix of the measurement model the same measurement.

    At each step, a state is proposed ``\\rho_{\\text{prop}} \\sim q(\\rho_{\\bullet,\\rho_{\\text{curr}})`` and the state is accepted with probability ``\\mathsf{P}(\\rho_i) = \\min(1,e^{-\\Delta_{i}})``

    input:

      - ``\\rho_{\\test{start}}``: initial state of the chain of states
      - ``\\Gamma_{\\text{sqrt}}``: square root matrix of the covariance of the smooth mechanism ``q_{\\text{smooth}}``
      - ``y``: measurement vector
      - ``y_d``: regularization expected values
      - ``\\Gamma_I^{-1}``: inverse covariance matrix of the measurement
      - ``\\Gamma_D^{-1}``: inverse covariance matrix of the regularization
      - ``H``: measurement operator
      - ``D``: regularization operator
      - ``N_s``: chain length
      - ``\\Gamma_{H,y}``: operator-to-noise covariance matrix

    output: 

      - ``(\\rho_i)_i``: chain of states
      - ``(\\Delta_{i})_{i}``: cost function variation associated to the chain of states
"""
function samplePosteriorMargin(ρ_start::Array{Cdouble,1},Γsqrt::Array{Cdouble,2},y::Array{Cdouble,1},yd::Array{Cdouble,1},ΓIinv::Array{Cdouble,2},Γdinv::Array{Cdouble,2},H::Array{Cdouble,2},D::Array{Cdouble,2},ΓHΓyinv::Array{Cdouble,2};Ns::Int64=10000)
    # all samples
    ρ_all = zeros(Cdouble,Ns+1,length(ρ_start))
    ρ_all[1,:] = ρ_start;
    deltaU = zeros(Cdouble,Ns);
    for i in 1:Ns
        # draw a new sample from a distribution not to far from the actual one
        ρ_all[i+1,:] = transmissionMechanismSmooth(ρ_all[i,:],Γsqrt)
        ρ_all[i+1,ρ_all[i+1,:].<0.0] .= 0.0  
        
        # accept or reject the sample
        ρ_all[i+1,:],deltaU[i] = acceptSampleMargin(ρ_all[i,:],ρ_all[i+1,:],y,yd,ΓIinv,Γdinv,H,D,ΓHΓyinv)
    end
    ρ_all,deltaU
end

