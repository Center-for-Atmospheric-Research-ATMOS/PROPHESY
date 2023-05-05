# push!(LOAD_PATH,"/home/mattoz/Dev/XPSinv.jl/src/")
# sudo nano /opt/julias/julia-1.3.1/etc/julia/startup.jl

"""

    This is the [`XPSinv`](@ref), which contains
    * [`XPSinv.D0th`](@ref)
    * [`XPSinv.D1st`](@ref)
    * [`XPSinv.D2nd`](@ref)
    * [`XPSinv.D3rd`](@ref)
    * [`XPSinv.D4th`](@ref)
    * [`XPSinv.reg_inv`](@ref)
    * [`XPSinv.pseudo_inv`](@ref)
    * [`XPSinv.alg2_cp_quad`](@ref) /
    * [`XPSinv.alg2_cp_quad_un`](@ref) /
    * [`XPSinv.alg2_cp_quad_LM`](@ref) /

    Derivation of a primal-dual optimization problem: implementation of ALG2, a.k.a. CP, in [1]
    [1] A. Chambolle and T. Pock, 2011. A first-order primal-dual algorithm for convex problems with applications to imaging.
    Journal of mathematical imaging and vision, 40(1), pp.120-145.
    DOI: 10.1007/s10851-010-0251-1

    Null Space optimization
    [2] D. Stolzenburg, M. Ozon, M.  Kulmala, K. Lehtinen, K. Lehtipalo and J. Kangasluoma, 2022, Combining instrument inversions for sub-10 nm aerosol number size-distribution measurements
    Journal of Aerosol Science, Vol. 159, p. 105862
    DOI: 10.1016/j.jaerosci.2021.105862

"""
module XPSinv

using Statistics
using LinearAlgebra #SVD and eigen value decomposition

# useful operators
export D0th,D1st, D2nd, D3rd, D4th, reg_inv, pseudo_inv


##
## ALG2
##

# call function for ALG2/CP
# export alg2_cp_gaussian, alg2_cp_gaussian_un, alg2_cp_gaussian_un_no_mem
# export alg2_cp_gaussian_un_val, alg2_cp_gaussian_un_no_mem_val
# export alg2_cp_gaussian_un_vals, alg2_cp_gaussian_un_no_mem_vals

export alg2_cp_quad, alg2_cp_quad_un # these are the two functions that lead to acceptable results (the only function of real use for this particular problem)
export alg2_cp_quad_LM               # not tracking all variable so that it limits the memory print

include("usualFun.jl")
# include("cp_gauss.jl")
include("cp_quad.jl")
include("cp_quad_uncertainty.jl")
# include("cp_gauss_uncertainty.jl")
# include("cp_gauss_uncertainty_bulk_value.jl")
# include("cp_gauss_uncertainty_values.jl")     # should replace cp_gauss_uncertainty_bulk_value.jl in future versions

# ##
# ## NSO
# ##
# export low_rank_reconstruction,rest_reconstruction,one_iteration_nso, null_nso, iterative_nso
# export low_rank_reconstruction_un, model_un_mat_img
# include("nso.jl")

# export low_rank_reconstruction_un_2, one_iteration_nso_un_2, null_nso_un_2, iterative_nso_un_2, data_sample_nso_un, model_un
# export shuffle_data, shuffle_data_sample_nso_un
# include("nso_un.jl")


end # module
