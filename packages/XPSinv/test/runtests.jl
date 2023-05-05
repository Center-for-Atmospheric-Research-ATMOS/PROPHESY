using Test
# using LinearAlgebra
using XPSinv



function test_D0th()
    # D0th
    true
end

function test_D1st()
    # D1st
    true
end

function test_D2nd()
    # D2nd
    true
end

function test_D3rd()
    # D3rd
    true
end

function test_D4th()
    # D4th
    true
end

function test_reg_inv()
    # reg_inv
    true
end

function test_pseudo_inv()
    # pseudo_inv
    true
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