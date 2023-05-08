""" 
0th order sparsity (meaning that a priori the signal should be 0 everywhere) 

input:

  - N: matrix dimension

output:

  - identity matrix
"""
function D0th(N::Int64)
    Matrix{Cdouble}(I,N,N)
end

""" 
first order sparsity: a priori constant 

```math
D_{i,j} = \\begin{cases}
1 & j=i\\\\
-1 & j=i+1\\\\
0 & \\text{otherwise}
\\end{cases}
```

input:

  - N: matrix dimension

output:

  - ``N-1\\times N`` matrix
"""
function D1st(N::Int64)
    diagm(N-1,N,0 => ones(Cdouble,N-1),  1 => -ones(Cdouble,N-1))
end

""" 
second order sparsity: a priori linear 

```math
D_{i,j} = \\begin{cases}
-1 & j=i,\\,j=i+2\\\\
2 & j=i+1\\\\
0 & \\text{otherwise}
\\end{cases}
```

input:

  - N: matrix dimension

output:

  - ``N-2\\times N`` matrix

"""
function D2nd(N::Int64)
    diagm(N-2,N,1 => 2ones(Cdouble,N-2), 0 => -ones(Cdouble,N-2) ,2 => -ones(Cdouble,N-2))
end

""" 
third order sparsity: a priori quadratic 

```math
D_{i,j} = \\begin{cases}
-1 & j=i\\\\
3 & j=i\\pm 1\\\\
-3 & j=i\\pm 2\\\\
1 & j=i\\pm 3\\\\
0 & \\text{otherwise}
\\end{cases}
```

input:

  - N: matrix dimension

output:

  - ``N-3\\times N`` matrix

"""
function D3rd(N::Int64)
    diagm(N-3,N,0 => -ones(Cdouble,N-3), 1 => 3.0ones(Cdouble,N-3) ,2 => -3.0ones(Cdouble,N-3) ,3 => ones(Cdouble,N-3))
end

""" 
forth order sparsity: a priori cubic 

```math
D_{i,j} = \\begin{cases}
-1 & j=i\\\\
4 & j=i\\pm 1\\\\
-6 & j=i\\pm 2\\\\
4 & j=i\\pm 3\\\\
-1 & j=i\\pm 4\\\\
0 & \\text{otherwise}
\\end{cases}
```

input:

  - N: matrix dimension

output:

  - ``N-4\\times N`` matrix

"""
function D4th(N::Int64)
    diagm(N-4,N,0 => -ones(Cdouble,N-4), 1 => 4.0ones(Cdouble,N-4) ,2 => -6.0ones(Cdouble,N-4) ,3 => 4.0ones(Cdouble,N-4) ,4 => -1.0ones(Cdouble,N-4))
end

# compactness: it's quite fair to assume that for z=0 the concentration of any species wuld be 0
function Dcompact(N::Int64,z_min::Cdouble,z_max::Cdouble,z0::Cdouble,σ0::Cdouble)
    z = collect(range(z_min,z_max,length=N));
    diagm(N,logistic.(z.-z0,1.0,0.0,1.0/σ0))
end


""" 
    reg_inv(H::Array{Cdouble,2},D::Array{Cdouble,2};τ::Cdouble=1.0)

    computes a regularized pseudo inverse of the matrix H based on the optimization problem
    ```math
    \\min_x ||Hx-y||^2 + \\tau||Dx||^2 
    ```

    whose solution is formally defined by

    ```math
    \\hat{x} =  (H^T H + \\tau D^T D)̂^{-1} H^T y 
    ```

    input:

      - H: measurement operator
      - D: regularization operator
      - ``\\tau``: regularization strength

    output:

      - pseudo inverse matrix ``(H^T H + \\tau D^T D)̂^{-1} H^T``

"""
function reg_inv(H::Array{Cdouble,2},D::Array{Cdouble,2};τ::Cdouble=1.0)
    inv(H'*H + τ*D'*D)*H'
end

""" 
pseudo inverse 

pseudo_inv(H::Array{Cdouble,2},n_singular::Int64;th_singular::Cdouble=0.0)

computes the an approximation of the matrix H using the singular value decomposition ``H = U\\text{diagm}(S)V^T``

```math
H_{\\text{approx}} = U \\text{diagm}(\\tilde{S})V^t ,\\,
\\tilde{S}_{i} = \\begin{cases}
\\max(S_i,th_{\\text{singular}}\\max(S)) & i\\leqslant n_{\\text{singular}}\\\\
0 & \\text{otherwise}
\\end{cases}
```

and 

```math
H_{\\text{approx}}^{\\dagger} = (U \\text{diagm}(\\tilde{S}^{-1})V^t)^T ,\\,
\\tilde{S}_{i}^{-1} = \\begin{cases}
\\frac{1}{\\max(S_i,th_{\\text{singular}}\\max(S))} & i\\leqslant n_{\\text{singular}}\\\\
0 & \\text{otherwise}
\\end{cases}
```

input:

  - H: measurement operator
  - n_singular: number of singular values kept in the spectral truncation
  - th_singular: relative threshold of the singular values

output:

  - approximation of the matrix H with n_singular non-zero singular values and a spectral threshold set to th_singular*maximum(singular value)
  - pseudo inverse with the same characteristic as the approximation
"""
function pseudo_inv(H::Array{Cdouble,2},n_singular::Int64;th_singular::Cdouble=0.0)
    # singular value decomposition
    F = svd(H, full=true);

    # bound away from 0
    Nke,N = size(H);
    W  = zeros(Cdouble,N,Nke);
    F.S[F.S.<(th_singular*maximum(F.S))] .= th_singular*maximum(F.S);
    W[1:n_singular,1:n_singular] = diagm(F.S[1:n_singular]); 
    Winv  = zeros(Cdouble,N,Nke);
    Winv[1:n_singular,1:n_singular] = diagm(1.0 ./F.S[1:n_singular]);

    # approximations of H and its pseudo inverse respectively
    (F.U*W')*F.Vt , ((F.U*Winv')*F.Vt)'
end
