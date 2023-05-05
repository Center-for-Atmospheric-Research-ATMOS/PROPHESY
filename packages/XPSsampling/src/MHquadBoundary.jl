"""
    acceptSampleBoundary(ρ_cur::Array{Cdouble,1},ρ_prop::Array{Cdouble,1},y::Array{Cdouble,1},ΓIinv::Array{Cdouble,2},H::Array{Cdouble,2},Dprior::Array{Cdouble,2},ρB::Array{Cdouble,1},σB::Array{Cdouble,1})

    Checks out the proposed state ``\\rho_{\\text{prop}}`` compared with the current state ``\\rho_{\\text{curr}}`` in terms of the energy function ``U`` to decide which state to accept.
    The proposed state is accepted with probability ``\\masthsf{P}(\\rho_{\\text{prop}}) = \\min(1,e^{-\\Delta})`` where ``\\Delta = U(\\rho_{\\text{prop}}) - U(\\rho_{\\text{curr}})``

    The cost function ``U`` is 
    ```math
        \\begin{equation}
            U(\\rho) = ||y-H\\rho||_{\\Gamma_I}^2 + ||y_d-D\\rho||_{\\Gamma_D}^2 + \\frac{1}{\\sigma_B^2} (\\rho^1-\\rho_B^1)^2 + \\frac{1}{\\sigma_B^2} (\\rho^N-\\rho_B^N)^2
        \\end{equation}
    ```
    where ``\\|y-H\\rho\\|_{\\Gamma_I}^2`` is the data fidelity term with ``y`` the vector of measurement data, ``H`` the measurement operator and ``\\Gamma_I`` the noise covariance matrix.
    The term ``\\|y_d-D\\rho\\|_{\\Gamma_D}^2`` is the regularizer of the inverse problem associated to the cost function ``U`` with ``D`` the regularization operator (e.g. second order difference), ``y_d`` a vector of expected values of the regularization and ``\\Gamma_D`` the covariance matrix of the regularization (strength and correlation)
    The last two terms force the estimate ``\\rho`` to the predefined values ``\\rho_B``.

    input:

      - ``\\rho_{\\test{cur}}`` and ``\\rho_{\\test{prop}}``: current and proposed states
      - ``y``: measurement vector
      - ``\\Gamma_I^{-1}``: inverse covariance matrix of the measurement
      - ``D^T\\Gamma_D^{-1}D``: product of the regularization opertor and inverse covariance matrix of the regularization
      - ``H``: measurement operator
      - ``D``: regularization operator
      - ``\\rho_B``: boundary values
      - ``\\sigma_B``: variances of the boundary values

    output: 

      - ``\\rho_{\\text{new}}``: the new state
      - ``-\\Delta``: the cost function variation
"""
function acceptSampleBoundary(ρ_cur::Array{Cdouble,1},ρ_prop::Array{Cdouble,1},y::Array{Cdouble,1},ΓIinv::Array{Cdouble,2},H::Array{Cdouble,2},Dprior::Array{Cdouble,2},ρB::Array{Cdouble,1},σB::Array{Cdouble,1})
    # if the posterior probability is larger for the proposed state ρ_prop than the current state ρ_cur, then accept the state, otherwise, reject it with probability p

    # likelihood
    r_cp = 0.5*(y-H*ρ_cur)'*ΓIinv*(y-H*ρ_cur)-0.5*(y-H*ρ_prop)'*ΓIinv*(y-H*ρ_prop);
    # smoothness
    r_cp = r_cp + 0.5ρ_cur'Dprior*ρ_cur - 0.5ρ_prop'Dprior*ρ_prop;
    # known values 
    r_cp = r_cp + 0.5*((ρ_cur[1]  -ρB[1])^2)/(σB[1]^2) - 0.5*((ρ_prop[1]  -ρB[1])^2)/(σB[1]^2);
    r_cp = r_cp + 0.5*((ρ_cur[end]-ρB[2])^2)/(σB[2]^2) - 0.5*((ρ_prop[end]-ρB[2])^2)/(σB[2]^2);
    
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
    acceptSampleModelMargin(ρ_cur::Array{Cdouble,1},ρ_prop::Array{Cdouble,1},y::Array{Cdouble,1},ΓIinv::Array{Cdouble,2},H::Array{Cdouble,2},ΓH::Array{Cdouble,3},Dprior::Array{Cdouble,2},ρB::Array{Cdouble,1},σB::Array{Cdouble,1})

    Checks out the proposed state ``\\rho_{\\text{prop}}`` compared with the current state ``\\rho_{\\text{curr}}`` in terms of the energy function ``U`` to decide which state to accept.
    The proposed state is accepted with probability ``\\masthsf{P}(\\rho_{\\text{prop}}) = \\min(1,e^{-\\Delta})`` where ``\\Delta = U(\\rho_{\\text{prop}}) - U(\\rho_{\\text{curr}})``

    The cost function ``U`` is 
    ```math
        \\begin{equation}
            U(\\rho) = ||y-H\\rho||_{\\Gamma_I}^2 + ||y_d-D\\rho||_{\\Gamma_D}^2 + \\frac{1}{\\sigma_B^2} (\\rho^1-\\rho_B^1)^2 + \\frac{1}{\\sigma_B^2} (\\rho^N-\\rho_B^N)^2 + \\underset{k}{\\sum} \\frac{\\rho^T \\Gamma_{H_k}\\rho}{\\sigma_k^2}     
        \\end{equation}
    ```
    where ``\\|y-H\\rho\\|_{\\Gamma_I}^2`` is the data fidelity term with ``y`` the vector of measurement data, ``H`` the measurement operator and ``\\Gamma_I`` the noise covariance matrix.
    The term ``\\|y_d-D\\rho\\|_{\\Gamma_D}^2`` is the regularizer of the inverse problem associated to the cost function ``U`` with ``D`` the regularization operator (e.g. second order difference), ``y_d`` a vector of expected values of the regularization and ``\\Gamma_D`` the covariance matrix of the regularization (strength and correlation)
    The two terms ``\\frac{1}{\\sigma_B}(\\rho^i-\\rho_B^i)^2`` force the estimate ``\\rho`` to the predefined values ``\\rho_B``.
    The last term ``\\underset{k}{\\sum} \\frac{\\rho^T \\Gamma_{H_k}\\rho}{\\sigma_k^2}`` is the approximation of the marginalization of errors in the measurement model (the noise variance of the channel ``k`` is ``\\sigma_k^2`` and the covariance matrix of the corresponding channel is ``\\Gamma_{H_k}``)

    input:

      - ``\\rho_{\\test{cur}}`` and ``\\rho_{\\test{prop}}``: current and proposed states
      - ``y``: measurement vector
      - ``\\Gamma_I^{-1}``: inverse covariance matrix of the measurement
      - ``(\\Gamma_{H_k})_k``: covariance matrixes of the measurement operator (one for each channel)
      - ``D^T\\Gamma_D^{-1}D``: product of the regularization opertor and inverse covariance matrix of the regularization
      - ``H``: measurement operator
      - ``D``: regularization operator
      - ``\\rho_B``: boundary values
      - ``\\sigma_B``: variances of the boundary values

    output: 

      - ``\\rho_{\\text{new}}``: the new state
      - ``-\\Delta``: the cost function variation
"""
function acceptSampleModelMargin(ρ_cur::Array{Cdouble,1},ρ_prop::Array{Cdouble,1},y::Array{Cdouble,1},ΓIinv::Array{Cdouble,2},H::Array{Cdouble,2},ΓH::Array{Cdouble,3},Dprior::Array{Cdouble,2},ρB::Array{Cdouble,1},σB::Array{Cdouble,1})
    # if the posterior probability is larger for the proposed state ρ_prop than the current state ρ_cur, then accept the state, otherwise, reject it with probability p

    # likelihood
    r_cp = 0.5*(y-H*ρ_cur)'*ΓIinv*(y-H*ρ_cur)-0.5*(y-H*ρ_prop)'*ΓIinv*(y-H*ρ_prop);
    # smoothness
    r_cp = r_cp + 0.5ρ_cur'Dprior*ρ_cur - 0.5ρ_prop'Dprior*ρ_prop;
    # known values 
    r_cp = r_cp + 0.5*((ρ_cur[1]  -ρB[1])^2)/(σB[1]^2) - 0.5*((ρ_prop[1]  -ρB[1])^2)/(σB[1]^2);
    r_cp = r_cp + 0.5*((ρ_cur[end]-ρB[2])^2)/(σB[2]^2) - 0.5*((ρ_prop[end]-ρB[2])^2)/(σB[2]^2);
    # marginalization of the model
    for i in 1:length(y)
        r_cp = r_cp + 0.5ΓIinv[i,i]*(ρ_cur'*ΓH[i,:,:]*ρ_cur - ρ_prop'*ΓH[i,:,:]*ρ_prop)
    end

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
    samplePosteriorBoundary(ρ_start::Array{Cdouble,1},Γsqrt::Array{Cdouble,2},y::Array{Cdouble,1},ΓIinv::Array{Cdouble,2},H::Array{Cdouble,2},Dprior::Array{Cdouble,2},ρB::Array{Cdouble,1},σB::Array{Cdouble,1};Ns::Int64=10000,psmooth::Cdouble=0.99)

    computes a series of samples ``(\\rho_i)_i`` drawn from the generative communication mechanism ``q(\\rho_{\\text{prop}},\\rho_{\\text{curr}})`` and the cost function variations ``(\\Delta_{i})_{i}`` associated with each move.
    The cost function ``U`` is 
    ```math
        \\begin{equation}
            U(\\rho) = ||y-H\\rho||_{\\Gamma_I}^2 + ||y_d-D\\rho||_{\\Gamma_D}^2 + \\frac{1}{\\sigma_B^2} (\\rho^1-\\rho_B^1)^2 + \\frac{1}{\\sigma_B^2} (\\rho^N-\\rho_B^N)^2
        \\end{equation}
    ```
    where ``\\|y-H\\rho\\|_{\\Gamma_I}^2`` is the data fidelity term with ``y`` the vector of measurement data, ``H`` the measurement operator and ``\\Gamma_I`` the noise covariance matrix.
    The term ``\\|y_d-D\\rho\\|_{\\Gamma_D}^2`` is the regularizer of the inverse problem associated to the cost function ``U`` with ``D`` the regularization operator (e.g. second order difference), ``y_d`` a vector of expected values of the regularization and ``\\Gamma_D`` the covariance matrix of the regularization (strength and correlation)
    The last two terms force the estimate ``\\rho`` to the predefined values ``\\rho_B``.

    At each step, a state is proposed ``\\rho_{\\text{prop}} \\sim q(\\rho_{\\bullet,\\rho_{\\text{curr}})`` and the state is accepted with probability ``\\masthsf{P}(\\rho_i) = \\min(1,e^{-\\Delta_{i}})``

    input:

      - ``\\rho_{\\test{start}}``: initial state of the chain of states
      - ``\\Gamma_{\\text{sqrt}}``: square root matrix of the covariance of the smooth mechanism ``q_{\\text{smooth}}``
      - ``y``: measurement vector
      - ``\\Gamma_I^{-1}``: inverse covariance matrix of the measurement
      - ``D^T\\Gamma_D^{-1}D``: product of the regularization opertor and inverse covariance matrix of the regularization
      - ``H``: measurement operator
      - ``\\rho_B``: boundary values
      - ``\\sigma_B``: variances of the boundary values
      - ``N_s``: chain length
      - ``p_{\\text{smooth}}``: the probability of choosing the smooth mechanism 

    output: 

      - ``(\\rho_i)_i``: chain of states
"""
function samplePosteriorBoundary(ρ_start::Array{Cdouble,1},Γsqrt::Array{Cdouble,2},y::Array{Cdouble,1},ΓIinv::Array{Cdouble,2},H::Array{Cdouble,2},Dprior::Array{Cdouble,2},ρB::Array{Cdouble,1},σB::Array{Cdouble,1};Ns::Int64=10000,psmooth::Cdouble=0.99)
    # all samples
    ρ_all = zeros(Cdouble,Ns+1,length(ρ_start))
    ρ_all[1,:] = ρ_start;
    for i in 1:Ns
        # draw a new sample from a distribution not to far from the actual one
        ρ_all[i+1,:] = transmissionMechanism(ρ_all[i,:],Γsqrt,σB;psmooth=psmooth)
        
        # accept or reject the sample
        ρ_all[i+1,:],_ = acceptSampleBoundary(ρ_all[i,:],ρ_all[i+1,:],y,ΓIinv,H,Dprior,ρB,σB)
    end
    ρ_all
end



"""
samplePosteriorModelMargin(ρ_start::Array{Cdouble,1},Γsqrt::Array{Cdouble,2},y::Array{Cdouble,1},ΓIinv::Array{Cdouble,2},H::Array{Cdouble,2},ΓH::Array{Cdouble,3},Dprior::Array{Cdouble,2},ρB::Array{Cdouble,1},σB::Array{Cdouble,1};Ns::Int64=10000,psmooth::Cdouble=0.99)

    computes a series of samples ``(\\rho_i)_i`` drawn from the generative communication mechanism ``q(\\rho_{\\text{prop}},\\rho_{\\text{curr}})`` and the cost function variations ``(\\Delta_{i})_{i}`` associated with each move.
    The cost function ``U`` is 
    ```math
        \\begin{equation}
            U(\\rho) = ||y-H\\rho||_{\\Gamma_I}^2 + ||y_d-D\\rho||_{\\Gamma_D}^2 + \\frac{1}{\\sigma_B^2} (\\rho^1-\\rho_B^1)^2 + \\frac{1}{\\sigma_B^2} (\\rho^N-\\rho_B^N)^2 + \\underset{k}{\\sum} \\frac{\\rho^T \\Gamma_{H_k}\\rho}{\\sigma_k^2}     
        \\end{equation}
    ```
    where ``\\|y-H\\rho\\|_{\\Gamma_I}^2`` is the data fidelity term with ``y`` the vector of measurement data, ``H`` the measurement operator and ``\\Gamma_I`` the noise covariance matrix.
    The term ``\\|y_d-D\\rho\\|_{\\Gamma_D}^2`` is the regularizer of the inverse problem associated to the cost function ``U`` with ``D`` the regularization operator (e.g. second order difference), ``y_d`` a vector of expected values of the regularization and ``\\Gamma_D`` the covariance matrix of the regularization (strength and correlation)
    The two terms ``\\frac{1}{\\sigma_B}(\\rho^i-\\rho_B^i)^2`` force the estimate ``\\rho`` to the predefined values ``\\rho_B``.
    The last term ``\\underset{k}{\\sum} \\frac{\\rho^T \\Gamma_{H_k}\\rho}{\\sigma_k^2}`` is the approximation of the marginalization of errors in the measurement model (the noise variance of the channel ``k`` is ``\\sigma_k^2`` and the covariance matrix of the corresponding channel is ``\\Gamma_{H_k}``)
    
    At each step, a state is proposed ``\\rho_{\\text{prop}} \\sim q(\\rho_{\\bullet,\\rho_{\\text{curr}})`` and the state is accepted with probability ``\\masthsf{P}(\\rho_i) = \\min(1,e^{-\\Delta_{i}})``

    input:

      - ``\\rho_{\\test{start}}``: initial state of the chain of states
      - ``\\Gamma_{\\text{sqrt}}``: square root matrix of the covariance of the smooth mechanism ``q_{\\text{smooth}}``
      - ``y``: measurement vector
      - ``\\Gamma_I^{-1}``: inverse covariance matrix of the measurement
      - ``(\\Gamma_{H_k})_k``: covariance matrixes of the measurement operator (one for each channel)
      - ``D^T\\Gamma_D^{-1}D``: product of the regularization opertor and inverse covariance matrix of the regularization
      - ``H``: measurement operator
      - ``\\rho_B``: boundary values
      - ``\\sigma_B``: variances of the boundary values
      - ``N_s``: chain length
      - ``p_{\\text{smooth}}``: the probability of choosing the smooth mechanism 

    output: 

      - ``(\\rho_i)_i``: chain of states
"""
function samplePosteriorModelMargin(ρ_start::Array{Cdouble,1},Γsqrt::Array{Cdouble,2},y::Array{Cdouble,1},ΓIinv::Array{Cdouble,2},H::Array{Cdouble,2},ΓH::Array{Cdouble,3},Dprior::Array{Cdouble,2},ρB::Array{Cdouble,1},σB::Array{Cdouble,1};Ns::Int64=10000,psmooth::Cdouble=0.99)
    # all samples
    ρ_all = zeros(Cdouble,Ns+1,length(ρ_start))
    ρ_all[1,:] = ρ_start;
    for i in 1:Ns
        # draw a new sample from a distribution not to far from the actual one
        ρ_all[i+1,:] = transmissionMechanism(ρ_all[i,:],Γsqrt,σB;psmooth=psmooth)
        
        # accept or reject the sample
        ρ_all[i+1,:],_ = acceptSampleModelMargin(ρ_all[i,:],ρ_all[i+1,:],y,ΓIinv,H,ΓH,Dprior,ρB,σB)
    end
    ρ_all
end
