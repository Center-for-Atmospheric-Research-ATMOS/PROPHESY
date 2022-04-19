# push!(LOAD_PATH,"/home/mattoz/Dev/XPSinv.jl/src/")
# sudo nano /opt/julias/julia-1.3.1/etc/julia/startup.jl

module XPSinv

using Statistics
using LinearAlgebra #SVD
# using Interpolations # used in XPSmeas.jl
# using Printf
# using utilsFun

# export the functions for the end user
export D0th,D1st, D2nd, D3rd, D4th, reg_inv, pseudo_inv
export F_gaussian, G_gaussian, F_convex_conjugate_gaussian, prox_F_conj_gaussian, prox_G_gaussian, alg2_cp_gaussian
export F_gaussian_un, G_gaussian_un, F_convex_conjugate_gaussian_un, G_convex_conjugate_gaussian_un, prox_F_conj_gaussian_un, prox_G_gaussian_un, alg2_cp_gaussian_un, alg2_cp_gaussian_un_no_mem
export F_gaussian_un_val, G_gaussian_un_val, F_convex_conjugate_gaussian_un_val, G_convex_conjugate_gaussian_un_val, prox_F_conj_gaussian_un_val, prox_G_gaussian_un_val, alg2_cp_gaussian_un_val, alg2_cp_gaussian_un_no_mem_val
export F_gaussian_un_vals, G_gaussian_un_vals, F_convex_conjugate_gaussian_un_vals, G_convex_conjugate_gaussian_un_vals, prox_F_conj_gaussian_un_vals, prox_G_gaussian_un_vals, alg2_cp_gaussian_un_vals, alg2_cp_gaussian_un_no_mem_vals

include("usualFun.jl")
include("cp_gauss.jl")
include("cp_gauss_uncertainty.jl")
include("cp_gauss_uncertainty_bulk_value.jl")
include("cp_gauss_uncertainty_values.jl")     # should replace cp_gauss_uncertainty_bulk_value.jl in future versions


export low_rank_reconstruction,rest_reconstruction,one_iteration_nso, null_nso, iterative_nso
export low_rank_reconstruction_un, model_un_mat_img
include("nso.jl")

export low_rank_reconstruction_un_2, one_iteration_nso_un_2, null_nso_un_2, iterative_nso_un_2, data_sample_nso_un, model_un
export shuffle_data, shuffle_data_sample_nso_un
include("nso_un.jl")


end # module
