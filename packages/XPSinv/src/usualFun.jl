# 0th order sparsity (meaning that a priori the signal should be 0 everywhere)
function D0th(N::Int64)
    Matrix{Cdouble}(I,N,N)
end

# first order sparsity: a priori constant
function D1st(N::Int64)
    # diagm(0 => ones(Cdouble,N-1),  1 => -ones(Cdouble,N-2))[1:end-1,:]
    diagm(N-1,N,0 => ones(Cdouble,N-1),  1 => -ones(Cdouble,N-1))
end

# second order sparsity: a priori linear
function D2nd(N::Int64)
    # D_2nd = diagm(0 => 2ones(Cdouble,N-1), -1 => -ones(Cdouble,N-2) ,1 => -ones(Cdouble,N-2));
    # D_2nd[1,:] = D_2nd[2,:];
    # D_2nd[end,:] = D_2nd[end-1,:];
    diagm(N-2,N,1 => 2ones(Cdouble,N-2), 0 => -ones(Cdouble,N-2) ,2 => -ones(Cdouble,N-2))
end

# third order sparsity: a priori quadratic
function D3rd(N::Int64)
    diagm(N-3,N,0 => -ones(Cdouble,N-3), 1 => 3.0ones(Cdouble,N-3) ,2 => -3.0ones(Cdouble,N-3) ,3 => ones(Cdouble,N-3))
end

# forth order sparsity: a priori cubic
function D4th(N::Int64)
    diagm(N-4,N,0 => -ones(Cdouble,N-4), 1 => 4.0ones(Cdouble,N-4) ,2 => -6.0ones(Cdouble,N-4) ,3 => 4.0ones(Cdouble,N-4) ,4 => -1.0ones(Cdouble,N-4))
end

# compactness: it's quite fair to assume that for z=0 the concentration of any species wuld be 0
function Dcompact(N::Int64,z_min::Cdouble,z_max::Cdouble,z0::Cdouble,σ0::Cdouble)
    z = collect(range(z_min,z_max,length=N));
    diagm(N,logistic.(z.-z0,1.0,0.0,1.0/σ0))
end

# min_x ||Hx-y||^2 + τ||Dx||^2 => x_opt =  (H^t H + τ D^t D)̂{-1} H^t y
function reg_inv(H::Array{Cdouble,2},D::Array{Cdouble,2};τ::Cdouble=1.0)
    inv(H'*H + τ*D'*D)*H'
end

function pseudo_inv(H::Array{Cdouble,2},n_singular::Int64;th_singular::Cdouble=0.0)
    # singular value decomposition
    F = svd(H, full=true);

    # bound away from 0
    Nke,N = size(H);
    W  = zeros(Cdouble,N,Nke);
    W[1:n_singular,1:n_singular] = diagm(th_singular*maximum(F.S) .+ F.S[1:n_singular]);
    Winv  = zeros(Cdouble,N,Nke);
    Winv[1:n_singular,1:n_singular] = diagm(1.0./(th_singular*maximum(F.S) .+ F.S[1:n_singular]));

    # approximations of H and its pseudo inverse respectively
    (F.U*W')*F.Vt , ((F.U*Winv')*F.Vt)'
end
