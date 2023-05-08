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
    # alg2_cp_quad
    true
end

function test_alg2_cp_quad_un()
    # alg2_cp_quad_un
    true
end

function test_alg2_cp_quad_LM()
    # alg2_cp_quad_LM
    true
end


@testset "Implementation of Pock and Chambolle algorithm" begin
    @test test_alg2_cp_quad()
    @test test_alg2_cp_quad_un()
    @test test_alg2_cp_quad_LM()
end