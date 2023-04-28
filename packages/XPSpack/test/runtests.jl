using Test
using XPSpack

function test1()
  true
end

function test2()
  true
end

@testset "XPSpack test bench" begin
  @test test1()
  @test test2()
end
