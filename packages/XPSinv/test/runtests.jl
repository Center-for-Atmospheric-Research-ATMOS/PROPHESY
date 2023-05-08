using Test
using LinearAlgebra
using XPSinv



function test_D0th()
    # D0th
    sum(D0th(10))==10.0
end

function test_D1st()
    # D1st
    sum(D1st(10))==0.0
end

function test_D2nd()
    # D2nd
    sum(D2nd(10))==0.0
end

function test_D3rd()
    # D3rd
    sum(D3rd(10))==0.0
end

function test_D4th()
    # D4th
    sum(D4th(10))==0.0
end

function test_reg_inv()
    N = 20;
    H = diagm(N-3,N,0 => ones(Cdouble,N-3),  1 => 0.5ones(Cdouble,N-3),  2 => 0.25ones(Cdouble,N-3),  3 => 0.125ones(Cdouble,N-3));
    D = D3rd(N);
    τ = 1.0;
    # reg_inv
    H_inv_reg = reg_inv(H,D;τ=τ)
    
    #results
    ((N,N-3)==size(H_inv_reg)) & (!isnan(sum(H_inv_reg))) & (!isinf(sum(H_inv_reg)))
end

function test_pseudo_inv()
    N = 25;
    H = diagm(N-3,N,0 => ones(Cdouble,N-3),  1 => 0.5ones(Cdouble,N-3),  2 => 0.25ones(Cdouble,N-3),  3 => 0.125ones(Cdouble,N-3));
    n_singular = 5;
    th_singular = 0.1;
    # pseudo_inv
    H_approx,H_inv = pseudo_inv(H,n_singular;th_singular=th_singular);
    
    # conditions
    cond1 = ((size(H)==size(H_approx)) & (!isnan(sum(H_approx))) & (!isinf(sum(H_approx))))
    cond2 = (((size(H,2),size(H,1))==size(H_inv)) & (!isnan(sum(H_inv))) & (!isinf(sum(H_inv))))

    # results
    cond1 & cond2
end

@testset "Utils functions" begin
    @test test_D0th()
    @test test_D1st()
    @test test_D2nd()
    @test test_D3rd()
    @test test_D4th()
    @test test_reg_inv()
    @test test_pseudo_inv()
end

function test_alg2_cp_quad()
    N = 20;
    # measurement operator
    H = diagm(N-3,N,0 => ones(Cdouble,N-3),  1 => 0.5ones(Cdouble,N-3),  2 => 0.25ones(Cdouble,N-3),  3 => 0.125ones(Cdouble,N-3));
    # regularization operator 
    D = D3rd(N);
    # augmented state operator
    A = [H;D];
    # data
    x = collect(LinRange(0.0,1.0,N));
    ρ = (x.-0.25).*(x.-0.75).+0.1;
    Γy = (0.01^2)*diagm(ones(Cdouble,N-3)); # measurement noise covariance matrix
    Γy_sqrt = sqrt(Γy);
    yy = H*ρ + Γy_sqrt*randn(N-3); # measurement data
    Γd = (0.01^2)*diagm(ones(Cdouble,N-3)); # covariance matrix of the a priori (strength and correlation)
    yd = zeros(Cdouble,N-3); # a priori data

    # algo setup
    ρ0 = zeros(Cdouble,N);    # initial conditions
    W_stop = ones(Cdouble,N); # weights for the stopping criteria
    r_n_tol = 1.0e-6;         # relative stopping criteria in the primal space
    r_y_tol = 0.5;            # relative stopping criteria in the data space
    Niter = 1000;             # maximum number of iteration
    τ0 = 1.0;                 # initial value of algorithm parameter τ (see Pock and Chambolle)

    # alg2_cp_quad
    xn,sn,τn,X_ALL,S_ALL,T_ALL,N_last = alg2_cp_quad(ρ0,yy,yd,A,Γy,Γd,W_stop;τ0=τ0,Niter=Niter,r_n_tol=r_n_tol,r_y_tol=r_y_tol);
    
    # conditions
    cond1 = (length(xn)==N) & (!isnan(sum(xn))) & (!isinf(sum(xn)))
    cond2 = (length(sn)==size(A,1)) & (!isnan(sum(sn))) & (!isinf(sum(sn)))
    cond3 = (!isnan(τn)) & (!isinf(τn))
    cond4 = ((N,Niter+1)==size(X_ALL)) & (!isnan(sum(X_ALL))) & (!isinf(sum(X_ALL)))
    cond5 = ((size(A,1),Niter+1)==size(S_ALL)) & (!isnan(sum(S_ALL))) & (!isinf(sum(S_ALL)))
    cond6 = (length(T_ALL)==(Niter+1)) & (!isnan(sum(T_ALL))) & (!isinf(sum(T_ALL)))
    cond7 = (N_last<=(Niter+1))

    # results
    cond1 & cond2 & cond3 & cond4 & cond5 & cond6 & cond7
end

function test_alg2_cp_quad_un()
    N = 20;
    # measurement operator
    H = diagm(N-3,N,0 => ones(Cdouble,N-3),  1 => 0.5ones(Cdouble,N-3),  2 => 0.25ones(Cdouble,N-3),  3 => 0.125ones(Cdouble,N-3));
    # global covariance of the uncertainty in the measurement operator
    ΓHΓyinv = 0.01diagm(ones(Cdouble,N));
    # regularization operator 
    D = D3rd(N);
    # augmented state operator
    A = [H;D;diagm(ones(Cdouble,N))];
    # data
    x = collect(LinRange(0.0,1.0,N));
    ρ = (x.-0.25).*(x.-0.75).+0.1;
    Γy = (0.01^2)*diagm(ones(Cdouble,N-3)); # measurement noise covariance matrix
    Γy_sqrt = sqrt(Γy);
    yy = H*ρ + Γy_sqrt*randn(N-3); # measurement data
    Γd = (0.01^2)*diagm(ones(Cdouble,N-3)); # covariance matrix of the a priori (strength and correlation)
    yd = zeros(Cdouble,N-3); # a priori data

    # algo setup
    ρ0 = zeros(Cdouble,N);    # initial conditions
    W_stop = ones(Cdouble,N); # weights for the stopping criteria
    r_n_tol = 1.0e-6;         # relative stopping criteria in the primal space
    r_y_tol = 0.5;            # relative stopping criteria in the data space
    Niter = 1000;             # maximum number of iteration
    τ0 = 1.0;                 # initial value of algorithm parameter τ (see Pock and Chambolle)

    # alg2_cp_quad_un
    xn,sn,τn,X_ALL,S_ALL,T_ALL,N_last = alg2_cp_quad_un(ρ0,yy,yd,A,Γy,Γd,ΓHΓyinv::Array{Cdouble,2},W_stop;τ0=τ0,Niter=Niter,r_n_tol=r_n_tol,r_y_tol=r_y_tol);
    
    # conditions
    cond1 = (length(xn)==N) & (!isnan(sum(xn))) & (!isinf(sum(xn)))
    cond2 = (length(sn)==size(A,1)) & (!isnan(sum(sn))) & (!isinf(sum(sn)))
    cond3 = (!isnan(τn)) & (!isinf(τn))
    cond4 = ((N,Niter+1)==size(X_ALL)) & (!isnan(sum(X_ALL))) & (!isinf(sum(X_ALL)))
    cond5 = ((size(A,1),Niter+1)==size(S_ALL)) & (!isnan(sum(S_ALL))) & (!isinf(sum(S_ALL)))
    cond6 = (length(T_ALL)==(Niter+1)) & (!isnan(sum(T_ALL))) & (!isinf(sum(T_ALL)))
    cond7 = (N_last<=(Niter+1))

    # results
    cond1 & cond2 & cond3 & cond4 & cond5 & cond6 & cond7
end

function test_alg2_cp_quad_LM()
    N = 20;
    # measurement operator
    H = diagm(N-3,N,0 => ones(Cdouble,N-3),  1 => 0.5ones(Cdouble,N-3),  2 => 0.25ones(Cdouble,N-3),  3 => 0.125ones(Cdouble,N-3));
    # regularization operator 
    D = D3rd(N);
    # augmented state operator
    A = [H;D];
    # data
    x = collect(LinRange(0.0,1.0,N));
    ρ = (x.-0.25).*(x.-0.75).+0.1;
    Γy = (0.01^2)*diagm(ones(Cdouble,N-3)); # measurement noise covariance matrix
    Γy_sqrt = sqrt(Γy);
    yy = H*ρ + Γy_sqrt*randn(N-3); # measurement data
    Γd = (0.01^2)*diagm(ones(Cdouble,N-3)); # covariance matrix of the a priori (strength and correlation)
    yd = zeros(Cdouble,N-3); # a priori data

    # algo setup
    ρ0 = zeros(Cdouble,N);    # initial conditions
    W_stop = ones(Cdouble,N); # weights for the stopping criteria
    r_n_tol = 1.0e-6;         # relative stopping criteria in the primal space
    r_y_tol = 0.5;            # relative stopping criteria in the data space
    Niter = 1000;             # maximum number of iteration
    τ0 = 1.0;                 # initial value of algorithm parameter τ (see Pock and Chambolle)

    # alg2_cp_quad_LM
    xn,sn,N_last = alg2_cp_quad_LM(ρ0,yy,yd,A,Γy,Γd,W_stop;τ0=τ0,Niter=Niter,r_n_tol=r_n_tol,r_y_tol=r_y_tol);
    
    # conditions
    cond1 = (length(xn)==N) & (!isnan(sum(xn))) & (!isinf(sum(xn)))
    cond2 = (length(sn)==size(A,1)) & (!isnan(sum(sn))) & (!isinf(sum(sn)))
    cond3 = (N_last<=(Niter+1))

    # results
    cond1 & cond2 & cond3
end


@testset "Implementation of Pock and Chambolle algorithm" begin
    @test test_alg2_cp_quad()
    @test test_alg2_cp_quad_un()
    @test test_alg2_cp_quad_LM()
end