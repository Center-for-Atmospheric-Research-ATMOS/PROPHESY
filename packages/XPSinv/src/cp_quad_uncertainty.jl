#------------------------------------------------------------------------------
#
# This file is part of the XPSinv module which is licensed under CC-BY 4.0 license.
#
# Copyright (C) 2022,  Matthew Ozon, Konstantin Tumashevich and Nønne L. Prisle.
#
#------------------------------------------------------------------------------

##
## proximal for gaussian noise problems
##

# z: argument of the cost function
# y: Gaussian distributed data
"""
    F_quad_un(z::Array{Cdouble,1},yy::Array{Cdouble,1},yd::Array{Cdouble,1},Γyinv::Array{Cdouble,2},Γdinv::Array{Cdouble,2},ΓHΓyinv::Array{Cdouble,2})

    Computes the value of the cost function for a quadratically-regularized gaussian cost
    function. The M first elements of the array z are the expectations of the
    Gaussian distribution related to measurement, the ``N_d`` next entries are the 
    expectation of the regularization and the last N entries are related to the marginalization
    of the model's error. The overall cost is quadratic and looks like:
        ```math
        F(z) = (yy-z_{1:M})^T \\Gamma_y^{-1} (yy-z_{1:M}) + (yd-z_{M+1:M+N_d})^T \\Gamma_d^{-1} (yd-z_{M+1:M+N_d}) + z_{M+N_d+1:M+N_d+N}^T \\Gamma_{H,y} z_{M+N_d+1:M+N_d+N}
        ```
    
    input:

      - z:       dual state
      - yy:      measurement data (length(yy)=``M``)
      - yd:      regularization data (length(yd)=``N_d``)
      - Γyinv:   inverse of the measurement data covariance matrix
      - Γdinv:   inverse of the covariance matrix of the regularization
      - ΓHΓyinv: operator-to-noise covariance matrix

    output:

      - evaluation of F(z)
"""
function F_quad_un(z::Array{Cdouble,1},yy::Array{Cdouble,1},yd::Array{Cdouble,1},Γyinv::Array{Cdouble,2},Γdinv::Array{Cdouble,2},ΓHΓyinv::Array{Cdouble,2})
    M = length(yy)
    N = length(yd)
    (yy-z[1:M])'*Γyinv*(yy-z[1:M]) + (yd-z[M+1:M+N])'*Γdinv*(yd-z[M+1:M+N]) + z[M+N+1:end]'*ΓHΓyinv*z[M+N+1:end]  # separate into data space, regularization space and model uncertainty space
end

"""
    G_quad_un(x::Array{Cdouble,1})

    computes the indicator function of the positive quadrant ``(\\mathbb{R}^+)^N``. It returns
    0 if x is in the positive quadrant and ``\\infty`` otherwise

    input:

        - x: primal state vector
"""
function G_quad_un(x::Array{Cdouble,1})
    val = 0.0
    if any(x.<0.0)
        val = Inf
    end
    val
end


## convex conjugate of F (convex conjugate = Legendre–Fenchel transformation = generalization of the Legendre transformation)
# F⋆(u) = sup_{z} ∑_i u[i]*z[i] - F(z)
# ∑_{i=1}^{M+N} u[i]y[i] + (1/4)*u[1:M]'*Γy*u[1:M] + (1/4)*u[M+1:M+N]'*Γd*u[M+1:M+N] + (1/4)*u[M+N+1:end]'*ΓHΓyinv^{-1}*u[M+N+1:end]
# u: argument of the convex conjugate of the cost function
# y: Gaussian distributed data
"""
    F_convex_conjugate_quad_un(u::Array{Cdouble,1},yy::Array{Cdouble,1},yd::Array{Cdouble,1},Γy::Array{Cdouble,2},Γd::Array{Cdouble,2},ΓHΓyinv::Array{Cdouble,2})

    Computes the convex conjugate of
        ```math
        F(z) = (yy-z_{1:M})^T \\Gamma_y^{-1} (yy-z_{1:M}) + (yd-z_{M+1:M+N_d})^T \\Gamma_d^{-1} (yd-z_{M+1:M+N_d}) + z_{M+N_d+1:M+N_d+N}^T \\Gamma_{H,y} z_{M+N_d+1:M+N_d+N}
        ```
    defined by:
        ```math
        F_{\\star}(u) = \\sup_z u^Tz - F(z)
        ```
    which, in the quadratically-regularized Gaussian case is given by:
        ```math
        \\sum_{i=1}^{M+N_d+N} u_i y_i + \\frac{1}{4} u_{1:M}^T \\Gamma_y u_{1:M} + \\frac{1}{4} u_{M+1:M+N_d}^T \\Gamma_d u_{M+1:M+N_d} + \\frac{1}{4} u_{M+N_d+1:M+N_d+N}^T \\Gamma_{H,y}^{-1} u_{M+N_D+1:M+N_d+N}
        ```
    where N is the dimension of the primal space

    input:

        - u:       dual state
        - yy:      measurement data (length(yy)=``M``)
        - yd:      regularization data (length(yd)=``N_d``)
        - Γy:      measurement data covariance matrix
        - Γd:      covariance matrix of the regularization
        - ΓHΓyinv: operator-to-noise covariance matrix
  
      output:
  
        - evaluation of the convex conjugate of F⋆ in u
"""
function F_convex_conjugate_quad_un(u::Array{Cdouble,1},yy::Array{Cdouble,1},yd::Array{Cdouble,1},Γy::Array{Cdouble,2},Γd::Array{Cdouble,2},ΓHΓyinv::Array{Cdouble,2})
    M = length(yy)
    N = length(yd)
    u[1:M]'*(yy + 0.25*Γy*u[1:M]) + u[M+1:M+N]'*(yd + 0.25*Γd*u[M+1:M+N]) + 0.25*u[M+N+1:end]'*inv(ΓHΓyinv)*u[M+N+1:end]
end

## proximal of the convex conjugate of F
"""
    prox_F_conj_quad_un(x::Array{Cdouble},yy::Array{Cdouble,1},yd::Array{Cdouble,1},σ::Cdouble,Py::Array{Cdouble,2},Dy::Array{Cdouble,1},Pd::Array{Cdouble,2},Dd::Array{Cdouble,1},Ph::Array{Cdouble,2},Dh::Array{Cdouble,1})

    Γy = Py*diagm(Dy)*Py'
    Γd = Pd*diagm(Dd)*Pd'
    ΓHΓyinv = Ph*diagm(Dh)*Ph'

    computes the proximal of σF⋆ evaluated in x, with F being defined by:
        ```math
        F(z) = (yy-z_{1:M})^T \\Gamma_y^{-1} (yy-z_{1:M}) + (yd-z_{M+1:M+N_d})^T \\Gamma_d^{-1} (yd-z_{M+1:M+N_d}) + z_{M+N_d+1:M+N_d+N}^T \\Gamma_{H,y} z_{M+N_d+1:M+N_d+N}
        ```
    and the proximal by:
        ```math
        \\text{prox}_{\\sigma F_{\\star}}(x) = \\arg\\min_u {\\sigma F_{\\star}(u) + \\frac{1}{2}*||u-x||^2}
        ```
    where N is the dimension of the primal space, M (length(yy)) the dimension of the measurment data and N_d (length(yd)) the dimension of the regularization data

    input:

    - x:       dual state (length(x)=``M+N_d+N``)
    - yy:      measurement data (length(yy)=``M``)
    - yd:      regularization data (length(yd)=``N_d``)
    - ``\\sigma``: proximal parameter
    - Py,Dy:   eigen decomposition matrices of the measurement data covariance matrix (``\\Gamma_y`` = Py*diagm(Dy)*Py')
    - Pd,Dd:   eigen decomposition matrices of the covariance matrix of the regularization (``\\Gamma_d`` = Pd*diagm(Dd)*Pd')
    - Ph,Dh:   eigen decomposition matrices of the operator-to-noise covariance matrix (``\\Gamma_{H,y}`` = Ph*diagm(Dh)*Ph')

    output:

      - proximal of the convex conjugate of the cost function F

"""
function prox_F_conj_quad_un(x::Array{Cdouble},yy::Array{Cdouble,1},yd::Array{Cdouble,1},σ::Cdouble,Py::Array{Cdouble,2},Dy::Array{Cdouble,1},Pd::Array{Cdouble,2},Dd::Array{Cdouble,1},Ph::Array{Cdouble,2},Dh::Array{Cdouble,1})
    M = length(yy)
    N = length(yd)
    u_data = Py*diagm(1.0 ./(1.0 .+ 0.5σ*Dy))*Py'*(x[1:M]-σ*yy)
    u_reg = Pd*diagm(1.0 ./(1.0 .+ 0.5σ*Dd))*Pd'*(x[M+1:M+N]-σ*yd)
    u_mod = Ph*diagm(1.0 ./(1.0 .+ 0.5σ./Dh))*Ph'*x[M+N+1:end]
    [u_data; u_reg; u_mod]
end

# proximal operator of the indicator function prox_{G}(x) = argmin_u {G(u) + (1/2)*||u-x||^2}
"""
    computes the proximal of G, the positive quadrant indicator function,
        ```math
        \\text{prox}_G(x) = \\arg\\min_u {G(u) + \\frac{1}{2}*||u-x||^2}
        ```
    input:

        - x: primal state vector
"""
function prox_G_quad_un(x::Array{Cdouble})
    x.*(x.>=0.0) # can't get any simpler
end

# x0:    initial point
# y:     data
# tau0:  initial value of τ in Alg2 chambolle2011first "A first-order primal-dual algorithm for convex problems with applications to imaging"
# Niter: number of iterations TODO: change the stopping criteria to a better one, maybe the relative variation between two steps or the distance between primal and dual
"""
    alg2_cp_quad_un(x0::Array{Cdouble,1},yy::Array{Cdouble,1},yd::Array{Cdouble,1},A::Array{Cdouble,2},Γy::Array{Cdouble,2},Γd::Array{Cdouble,2},ΓHΓyinv::Array{Cdouble,2},W_stop::Array{Cdouble,1};τ0::Cdouble=1.0,Niter::Int64=100,r_n_tol::Cdouble=1.0e-6,r_y_tol::Cdouble=0.5)

    run the algorithm 2 described in [1] for the minimization problem
        ```math
        \\min\\{F(Ax) + G(x)\\}
        ```
    with:
        ```math
        F(z) = (yy-z_{1:M})^T \\Gamma_y^{-1} (yy-z_{1:M}) + (yd-z_{M+1:M+N_d})^T \\Gamma_d^{-1} (yd-z_{M+1:M+N_d}) + z_{M+N_d+1:M+N_d+N}^T \\Gamma_{H,y} z_{M+N_d+1:M+N_d+N}
        ```
    and G the indicator function of the positive quadrant ``(\\mathbb{R}^+)^N``.
    N is the dimension of the primal space, M (length(yy)) the dimension of the measurment data and N_d (length(yd)) the dimension of the regularization data.

    input:

        - ``x_0``:              initial point (length(x0)=``N``)
        - ``y``:                measurement data vector (length(yy)=``M``)
        - ``y_d``:              vector of expected value of the regularization (length(yd)=``N_d``)
        - ``A``:                augmented measurement operator (measurement operator H, regularization operator D, and identity I, A=[H;D;I])
        - ``\\Gamma_y``:        covariance matrix of the measurement data
        - ``\\Gamma_d``:        covariance matrix of the regularization (strength and correlation)
        - ``\\Gamma_{H,y}``:    operator-to-noise covariance matrix
        - ``W_{\\text{stop}}``: weights used for the stopping criteria in the primal space
      
      optional input:
  
        - ``\\tau_0``:          initial value of algorithm parameter τ (see Pock and Chambolle [1]), default value ``\\tau_0=1``
        - ``N_{\\text{iter}}``: maximum number of iteration, default value ``N_{\\text{iter}}=100``
        - r_n_tol:              relative stopping criteria in the primal space, default value r_n_tol=1.0e-6
        - r_y_tol:              relative stopping criteria in the data space, default value r_y_tol=0.5
  
      output:
  
        - xn:                   final primal state 
        - sn:                   final dual state
        - τn:                   final value of the parameter ``\\tau``
        - X_ALL:                all primal states along the iterations
        - S_ALL:                all dual states along the iterations
        - T_ALL:                all values of ``\\tau`` along the iterations
        - N_last:               number of iteration computed (the algorithm stops either because it reached the maximum number of iteration or because a stopping criteria is met)
  

    [1] A. Chambolle and T. Pock, 2011. A first-order primal-dual algorithm for convex problems with applications to imaging.
    Journal of mathematical imaging and vision, 40(1), pp.120-145.
"""
function alg2_cp_quad_un(x0::Array{Cdouble,1},yy::Array{Cdouble,1},yd::Array{Cdouble,1},A::Array{Cdouble,2},Γy::Array{Cdouble,2},Γd::Array{Cdouble,2},ΓHΓyinv::Array{Cdouble,2},W_stop::Array{Cdouble,1};τ0::Cdouble=1.0,Niter::Int64=100,r_n_tol::Cdouble=1.0e-6,r_y_tol::Cdouble=0.5)
    # init algorithm's parameters
    N = length(x0);
    if (N!=size(A,2))
        throw("alg2_cp_quad_un: operator A and state x0 do not have compatible dimensions")
    end
    if ((N+length(yy)+length(yd))!=size(A,1))
        throw("alg2_cp_quad_un: operator A does not have the right number of rows")
    end
    if (N!=length(W_stop))
        throw("alg2_cp_quad_un: stopping criteria dimension")
    end
    if ((length(yy)!=size(Γy,1)) | (length(yy)!=size(Γy,2)))
        throw("alg2_cp_quad_un: measurement data and covariance matrix dimension")
    end
    if ((length(yd)!=size(Γd,1)) | (length(yd)!=size(Γd,2)))
        throw("alg2_cp_quad_un: regularization data and covariance matrix dimension")
    end
    if ((N!=size(ΓHΓyinv,1)) | (N!=size(ΓHΓyinv,2)))
        throw("alg2_cp_quad_un: operator-to-noise covariance matrix dimension")
    end
    L = opnorm(A,2); # spectral norm of the operator A (if it's a matrix, it is the Lipschitz constant)
    σ0 = 1.0/(τ0*L^2);
    dxhx0 = sqrt(N); # just a rough estimate (upper estimate) of the distance between the initial state and the true state ||x^{\star}-x0||
    γ = 2.0*dxhx0/τ0; # makes the convergence faster

    # eigen decomposition of the covariance matrices
    Fy = eigen(Γy)
    Py = real(Fy.vectors)
    Dy = real(Fy.values)
    Fd = eigen(Γd)
    Pd = real(Fd.vectors)
    Dd = real(Fd.values)
    Fh = eigen(ΓHΓyinv)
    Ph = real(Fh.vectors)
    Dh = real(Fh.values)
 
    X_ALL = zeros(Cdouble,length(x0),Niter+1)
    S_ALL = zeros(Cdouble,size(A,1),Niter+1)
    T_ALL = zeros(Cdouble,Niter+1)

    # iteration
    σn     = σ0;
    τn     = τ0;
    xn     = copy(x0);
    sn     = zeros(Cdouble,size(A,1)); # could this be a problem? is sn allowed to take whatever value?
    xn_pre = copy(x0);
    xn_bar = copy(x0);
    X_ALL[:,1] = xn
    S_ALL[:,1] = sn
    T_ALL[1]   = τn
    idxdd = findall(yy.>0)
    idxnn = findall(yy.==0)
    Y_rel = zeros(Cdouble,length(yy));
    dX_W = copy(x0);
    N_last = Niter
    for i in 1:Niter
        # dual step
        sn = prox_F_conj_quad_un(sn+σn*A*xn_bar,yy,yd,σn,Py,Dy,Pd,Dd,Ph,Dh)
        # primal step
        xn_pre = xn; 
        xn     = prox_G_quad_un(xn-τn*A'*sn);
        # update algortihm paramters
        θn = 1.0/sqrt(1.0+2.0γ*τn);
        τn   = θn*τn;
        σn = σn/θn;
        # bar update: linear combination of the previous and current state
        xn_bar = xn + θn*(xn-xn_pre);

        # compute stopping criteria
        #    data fidelity
        Y_rel[idxdd] = abs.(yy[idxdd] - A[idxdd,:]*xn)./yy[idxdd]
        Y_rel[idxnn] = abs.(A[idxnn,:]*xn)
        #    norm of the relative step length
        dX_W = W_stop.*(xn-xn_pre);
        normdX_W = sum(dX_W.^2);
        normX_W = sum((W_stop.*xn).^2)
        #    both criteria must be met: data fidelity and relative change
        if ((median(Y_rel)<=r_y_tol) & (sqrt.(normdX_W./normX_W)<=r_n_tol))
            N_last = i;
            i = Niter+1;
            break
        end

        X_ALL[:,i+1] = xn
        S_ALL[:,i+1] = sn
        T_ALL[i+1]   = τn
    end

    # return primal and dual states (maybe the primal and dual problems' values and the path of x and s for debugging purposes)
    xn,sn,τn,X_ALL,S_ALL,T_ALL,N_last
end
