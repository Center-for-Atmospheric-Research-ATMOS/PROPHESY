# push!(LOAD_PATH,"/home/mattoz/Dev/XPSinv.jl/src/")
# sudo nano /opt/julias/julia-1.3.1/etc/julia/startup.jl

# Null Space optimization
#     [2] D. Stolzenburg, M. Ozon, M.  Kulmala, K. Lehtinen, K. Lehtipalo and J. Kangasluoma, 2022, Combining instrument inversions for sub-10 nm aerosol number size-distribution measurements
#     Journal of Aerosol Science, Vol. 159, p. 105862
#     DOI: 10.1016/j.jaerosci.2021.105862

"""
This is the [`XPSinv`](@ref), which contains
* [`XPSinv.D0th`](@ref)
* [`XPSinv.D1st`](@ref)
* [`XPSinv.D2nd`](@ref)
* [`XPSinv.D3rd`](@ref)
* [`XPSinv.D4th`](@ref)
* [`XPSinv.reg_inv`](@ref)
* [`XPSinv.pseudo_inv`](@ref)
* [`XPSinv.alg2_cp_quad`](@ref)
* [`XPSinv.alg2_cp_quad_un`](@ref)
* [`XPSinv.alg2_cp_quad_LM`](@ref)
* [`XPSinv.F_quad`](@ref)
* [`XPSinv.G_quad`](@ref)
* [`XPSinv.F_convex_conjugate_quad`](@ref)
* [`XPSinv.prox_F_conj_quad`](@ref)
* [`XPSinv.prox_G_quad`](@ref)
* [`XPSinv.F_quad_un`](@ref)
* [`XPSinv.G_quad_un`](@ref)
* [`XPSinv.F_convex_conjugate_quad_un`](@ref)
* [`XPSinv.prox_F_conj_quad_un`](@ref)
* [`XPSinv.prox_G_quad_un`](@ref)

    Derivation of a primal-dual optimization problem: implementation of ALG2, a.k.a. CP, in [1]
    [1] A. Chambolle and T. Pock, 2011. A first-order primal-dual algorithm for convex problems with applications to imaging.
    Journal of mathematical imaging and vision, 40(1), pp.120-145.
    DOI: 10.1007/s10851-010-0251-1

"""
module XPSinv

using Statistics
using LinearAlgebra #SVD and eigen value decomposition

# useful operators
export D0th,D1st, D2nd, D3rd, D4th, reg_inv, pseudo_inv


# call function for ALG2/CP
export alg2_cp_quad, alg2_cp_quad_un # these are the two functions that lead to acceptable results (the only function of real use for this particular problem)
export alg2_cp_quad_LM               # not tracking all variable so that it limits the memory print

include("usualFun.jl")
include("cp_quad.jl")
include("cp_quad_uncertainty.jl")


end # module
