##
## proximal for gaussian noise problems (with some uncertainty approximation)
##

# z: argument of the cost function
# y: Gaussian distributed data
"""
    F_gaussian(z::Array{Cdouble,1},y::Array{Cdouble,1})

    Computes the value of the cost function for a quadratically-regularized gaussian cost
    function. The M(=length(y)) first elements of the array z are the expectations of the
    Gaussian distribution. The rest of the array z contains the regularization entries,
    e.g. z[M+1:end] = R*n, to be quadratically penalized. The overall cost is:
        F(z) = ∑_{i=1}^M (z[i] - y[i])^2 + ∑_{i=M+1}^{end} z[i]^2
"""
function F_gaussian_un_val(z::Array{Cdouble,1},y::Array{Cdouble,1},ρ0::Cdouble,ρB::Cdouble,γ::Array{Cdouble,1},γ_H::Array{Cdouble,2},γ_D::Array{Cdouble,1},γ0::Cdouble,γB::Cdouble)
    Nke = length(y);
    M = size(γ_H,2);
    M_d = length(γ_D);
    W_i = (1.0/4.0)./dropdims(sum(((1.0./γ).*γ_H).^2,dims=1),dims=1) # should be of dimension M
    sum(((z[1:Nke]-y)./γ).^2) + (1.0./γ)'*((γ_H.^2)*(z[Nke+1:Nke+M].^2)) + sum((1.0./γ_D.^2).*z[Nke+M+1:Nke+M+M_d].^2) + ((z[end-1]-ρ0)/γ0)^2 + ((z[end]-ρB)/γB)^2
end

"""
    G_gaussian(x::Array{Cdouble,1})

    computes the indicator function of the positive quadrant (R^+)^N. It returns
    0 if x is in the positive quadrant and +∞ otherwise
"""
function G_gaussian_un_val(x::Array{Cdouble,1})
    val = 0.0
    if any(x.<0.0)
        val = Inf
    end
    val
end

## convex conjugate of F (convex conjugate = Legendre–Fenchel transformation = generalization of the Legendre transformation)
# u: argument of the convex conjugate of the cost function
# y: Gaussian distributed data
"""
    F_convex_conjugate_gaussian(u::Array{Cdouble,1},y::Array{Cdouble,1})

    Computes the convex conjugate of
        F(z) = ∑_{i=1}^M (z[i] - y[i])^2 + ∑_{i=M+1}^{end} z[i]^2
    defined by:
        F⋆(u) = sup_{z} ∑_i u[i]*z[i] - F(z)
    which, in the quadratically-regularized Gaussian case is given by:
        ∑_{i=1}^M u[i]y[i] + (1/4)*∑_{i=1}^{end} u[i]^2
"""
function F_convex_conjugate_gaussian_un_val(u::Array{Cdouble,1},y::Array{Cdouble,1},ρ0::Cdouble,ρB::Cdouble,γ::Array{Cdouble,1},γ_H::Array{Cdouble,2},γ_D::Array{Cdouble,1},γ0::Cdouble,γB::Cdouble)
    Nke = length(y);
    M = size(γ_H,2);
    M_d = length(γ_D);
    W_i = (1.0/4.0)./dropdims(sum(((1.0./γ).*γ_H).^2,dims=1),dims=1) # should be of dimension M
    sum(u[1:Nke].*(y+0.25(γ.^2).*u[1:Nke])) + sum(W_i.*u[Nke+1:Nke+M]) + sum(0.25(γ_D.^2).*u[Nke+M+1:Nke+M+M_d]) + u[end-1].*(ρ0+0.25(γ0.^2).*u[end-1]) + u[end].*(ρB+0.25(γB.^2).*u[end])
end

function G_convex_conjugate_gaussian_un_val(u::Array{Cdouble,1})
    val = 0.0
    if any(u.>0.0)
        val = Inf
    end
    val
end

## proximal of the convex conjugate of F
"""
    prox_F_conj(x::Array{Cdouble},y::Array{Cdouble,1},sigma::Cdouble)

    computes the proximal of σF⋆ evaluated in x, with F being defined by:
        F(z) = ∑_{i=1}^M (z[i] - y[i])^2 + ∑_{i=M+1}^{end} z[i]^2
    and the proximal by:
        prox_{σF⋆}(x) = argmin_u {σF⋆(u) + (1/2)*||u-x||^2}
"""
function prox_F_conj_gaussian_un_val(x::Array{Cdouble},y::Array{Cdouble,1},ρ0::Cdouble,ρB::Cdouble, σ::Cdouble,γ::Array{Cdouble,1},γ_H::Array{Cdouble,2},γ_D::Array{Cdouble,1},γ0::Cdouble,γB::Cdouble)
    Nke = length(y);
    M = size(γ_H,2);
    M_d = length(γ_D);

    W_i = (1.0/4.0)./dropdims(sum(((1.0./γ).*γ_H).^2,dims=1),dims=1)
    u_all = zeros(Cdouble,Nke+M+M_d+2);

    # in the data space
    for i in 1:Nke
        u_all[i] = (x[i]-σ*y[i])/(1.0+0.5σ*(γ[i]^2)); #
    end

    # in the model uncertainty space
    for i in 1:M
        u_all[Nke+i] = x[Nke+i]/(1.0+0.5σ*W_i[i])
    end

    # in the regularization space
    for i in 1:M_d
        u_all[Nke+M+i] = x[Nke+M+i]/(1.0+0.5σ*(γ_D[i]^2))
    end

    # the a priori known values
    u_all[end-1] = (x[end-1]-σ*ρ0)/(1.0+0.5σ*(γ0^2));
    u_all[end]   = (x[end]-σ*ρB)/(1.0+0.5σ*(γB^2));

    u_all
end

# proximal operator of the indicator function
"""
    computes the proximal of G, the positive quadrant indicator function,
        prox_{G}(x) = argmin_u {G(u) + (1/2)*||u-x||^2}
"""
function prox_G_gaussian_un_val(x::Array{Cdouble})
    x.*(x.>=0.0) # can't get any simpler
end

# x0:    initial point
# y:     data
# tau0:  initial value of \tau in Alg2 chambolle2011first "A first-order primal-dual algorithm for convex problems with applications to imaging"
# Niter: number of iterations TODO: change the stopping criteria to a better one, maybe the relative variation between two steps or the distance between primal and dual
"""
    alg2_cp_gaussian(x0::Array{Cdouble,1},y::Array{Cdouble,1},A::Array{Cdouble,2};tau0::Cdouble=1.0,Niter::Int64=100)

    run the algorithm 2 described in [1] for the minimization problem
        min {F(A*x) + G(x)}
    with:
        F(z) = ∑_{i=1}^M (z[i] - y[i])^2 + ∑_{i=M+1}^{end} z[i]^2
    and G the indicator function of the positive quadrant (R^+)^N.

    [1] A. Chambolle and T. Pock, 2011. A first-order primal-dual algorithm for convex problems with applications to imaging.
    Journal of mathematical imaging and vision, 40(1), pp.120-145.
"""
function alg2_cp_gaussian_un_val(x0::Array{Cdouble,1},s0::Array{Cdouble,1},y::Array{Cdouble,1},ρ0::Cdouble,ρB::Cdouble,A::Array{Cdouble,2},γ::Array{Cdouble,1},γ_H::Array{Cdouble,2},γ_D::Array{Cdouble,1},γ0::Cdouble,γB::Cdouble,W_stop::Array{Cdouble,1};tau0::Cdouble=1.0,Niter::Int64=100,r_n_tol::Cdouble=1.0e-6,r_y_tol::Cdouble=0.5)
    # init algorithm's parameters
    N = length(x0);
    if (N!=size(A,2))
        throw("alg2_chpo: operator A and state x0 do not have compatible dimensions")
    end
    L = opnorm(A,2); # spectral norm of the operator A (if it's a matrix, it is the Lipschitz constant)
    sigma0 = 1.0/(tau0*L^2);
    dxhx0 = sqrt(N); # sqrt(N)*1.4e4; # it is just a rough estimate (upper estimate) of the distance between the initial state and the true size distribution ||x^{\star}-x0||
    # gamma = 100.0*dxhx0/tau0;
    gamma = 2.0*dxhx0/tau0; # makes the convergence faster

    X_ALL = zeros(Cdouble,length(x0),Niter+1)
    S_ALL = zeros(Cdouble,size(A,1),Niter+1)
    T_ALL = zeros(Cdouble,Niter+1)

    # iteration
    sigman = sigma0;
    taun   = tau0
    xn     = copy(x0);
    sn     = copy(s0); # zeros(Cdouble,size(A,1)); # could this be a problem? is sn allowed to take whatever value?
    xn_pre = copy(x0);
    xn_bar = copy(x0);
    X_ALL[:,1] = xn
    S_ALL[:,1] = sn
    T_ALL[1]   = taun
    idxdd = findall(y.>0)
    idxnn = findall(y.==0)
    Y_rel = zeros(Cdouble,length(y));
    dX_W = copy(x0);
    N_last = Niter
    for i in 1:Niter
        # dual step
        sn = prox_F_conj_gaussian_un_val(sn+sigman*A*xn_bar,y,ρ0,ρB,sigman,γ,γ_H,γ_D,γ0,γB);
        # primal step
        xn_pre = xn; # memory (TODO: check that it does what I think it does)
        xn     = prox_G_gaussian_un_val(xn-taun*A'*sn);
        # update algortihm paramters
        thetan = 1.0/sqrt(1.0+2.0gamma*taun);
        taun   = thetan*taun;
        sigman = sigman/thetan;
        # bar update: linear combination of the previous and current state
        xn_bar = xn + thetan*(xn-xn_pre);

        # compute stopping criteria
        #    data fidelity
        Y_rel[idxdd] = abs.(y[idxdd] - A[idxdd,:]*xn)./y[idxdd]
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
        T_ALL[i+1]   = taun
    end

    # return primal and dual states (maybe the primal and dual problems' values and the path of x and s for debugging purposes)
    xn,sn,taun,X_ALL,S_ALL,T_ALL,N_last
end


function alg2_cp_gaussian_un_no_mem_val(x0::Array{Cdouble,1},s0::Array{Cdouble,1},y::Array{Cdouble,1},ρ0::Cdouble,ρB::Cdouble,A::Array{Cdouble,2},γ::Array{Cdouble,1},γ_H::Array{Cdouble,2},γ_D::Array{Cdouble,1},γ0::Cdouble,γB::Cdouble,W_stop::Array{Cdouble,1};tau0::Cdouble=1.0,Niter::Int64=100,r_n_tol::Cdouble=1.0e-6,r_y_tol::Cdouble=0.5)
    # init algorithm's parameters
    N = length(x0);
    if (N!=size(A,2))
        throw("alg2_chpo: operator A and state x0 do not have compatible dimensions")
    end
    L = opnorm(A,2); # spectral norm of the operator A (if it's a matrix, it is the Lipschitz constant)
    sigma0 = 1.0/(tau0*L^2);
    dxhx0 = sqrt(N); # sqrt(N)*1.4e4; # it is just a rough estimate (upper estimate) of the distance between the initial state and the true size distribution ||x^{\star}-x0||
    # gamma = 100.0*dxhx0/tau0;
    gamma = 2.0*dxhx0/tau0; # makes the convergence faster

    T_ALL = zeros(Cdouble,Niter+1);
    F_ALL = zeros(Cdouble,Niter+1);
    G_ALL = zeros(Cdouble,Niter+1);
    F_STAR_ALL = zeros(Cdouble,Niter+1);
    G_STAR_ALL = zeros(Cdouble,Niter+1);
    INNER_PROD = zeros(Cdouble,Niter+1);

    # iteration
    sigman = sigma0;
    taun   = tau0
    xn     = copy(x0);
    sn     = copy(s0); # zeros(Cdouble,size(A,1)); # could this be a problem? is sn allowed to take whatever value?
    xn_pre = copy(x0);
    xn_bar = copy(x0);
    idxnn = findall(y.!=0.0); #yup... it's not the best way to deal with zeros
    y_nz = y[idxnn];
    A_nz = A[idxnn,:];
    Y_rel = abs.(y_nz - A_nz*xn)./y_nz;
    dX_W = copy(x0);
    F_ALL[1]      = F_gaussian_un_val(A*xn,y,ρ0,ρB,γ,γ_H,γ_D,γ0,γB)
    G_ALL[1]      = G_gaussian_un_val(xn);
    F_STAR_ALL[1] = F_convex_conjugate_gaussian_un_val(sn,y,ρ0,ρB,γ,γ_H,γ_D,γ0,γB);
    G_STAR_ALL[1] = G_convex_conjugate_gaussian_un_val(-A'*sn);
    INNER_PROD[1] = (A*xn)'*sn;
    T_ALL[1]   = taun
    N_last = Niter
    for i in 1:Niter
        # dual step
        sn = prox_F_conj_gaussian_un_val(sn+sigman*A*xn_bar,y,ρ0,ρB,sigman,γ,γ_H,γ_D,γ0,γB);
        # primal step
        xn_pre = xn; # memory
        xn     = prox_G_gaussian_un_val(xn-taun*A'*sn);
        # update algortihm paramters
        thetan = 1.0/sqrt(1.0+2.0gamma*taun);
        taun   = thetan*taun;
        sigman = sigman/thetan;
        # bar update: linear combination of the previous and current state
        xn_bar = xn + thetan*(xn-xn_pre);

        # compute stopping criteria
        #    data fidelity
        Y_rel = abs.(y_nz - A_nz*xn)./y_nz;
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

        # compute the values for the cost function
        F_ALL[i+1]      = F_gaussian_un_val(A*xn,y,ρ0,ρB,γ,γ_H,γ_D,γ0,γB)
        G_ALL[i+1]      = G_gaussian_un_val(xn);
        F_STAR_ALL[i+1] = F_convex_conjugate_gaussian_un_val(sn,y,ρ0,ρB,γ,γ_H,γ_D,γ0,γB);
        G_STAR_ALL[i+1] = G_convex_conjugate_gaussian_un_val(-A'*sn);
        INNER_PROD[i+1] = (A*xn)'*sn;
        T_ALL[i+1]      = taun
    end

    # return primal and dual states (maybe the primal and dual problems' values and the path of x and s for debugging purposes)
    xn,sn,taun,F_ALL,G_ALL,F_STAR_ALL,G_STAR_ALL,INNER_PROD,T_ALL,N_last
end
