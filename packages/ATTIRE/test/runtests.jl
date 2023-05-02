# ATTIRE test bench

using Test
using ATTIRE


#######################################
##        define test functions      ##
#######################################

function test_φi()
    T = 1.0;
    ΔKi = 1.0
    Ki = 1.0
    Ke = 1.0
    (φi(Ke,Ki,ΔKi,T)==(1.0/sqrt(2π)))
end

function test_Φi()
    Ki = collect(0.0:0.05:1.0)
    ΔKi = 1.0
    T = 1.0
    f_array = Φi(Ki,ΔKi,T)
    cond = true
    [ cond = (cond & (f_array[i](Ki[i])==(1.0/sqrt(2π)))) for i in eachindex(Ki)]
    (typeof(f_array)==Vector{Function}) & (length(Ki)==length(f_array)) & cond
end

function test_sourceSpread()
    hν = 1.0
    hνj = 1.0
    Δνj = 1.0
    Fνj = 1.0
    sourceSpread(hν,hνj,Δνj,Fνj)==(1.0/sqrt(2π))
end

function test_countElectrons()
    σ = 0.0
    S = 4.5
    Ns = 10000;
    noise_samples = Array{Cdouble,1}(undef,Ns)
    [noise_samples[i] = countElectrons(S,σ) for i in 1:Ns]
    μ = sum(noise_samples)/Ns
    σ2 = sum((noise_samples.-μ).^2)/Ns

    (abs(σ2/μ -1.0)<0.1) # this test may fail with a very small probability
end

function test_T_r4000()
    Ke = 0.0
    E_pass = 1.0
    (T_r4000(Ke,E_pass)==1.0)
end

function test_simulateSpectrum()
    Fνj = 1.0
    hνj = 1000.0
    Δνj = 1.0;
    Ki = collect(499.0:0.1:501.0)
    ΔKi = 0.2
    T = 1.0
    Kdisc = collect(498.4:0.05:501.6)
    Be0 = collect(499.0:0.1:501.0)
    σ_cs_fg = ones(Cdouble,length(Be0))
    σ_bg_vec = zeros(Cdouble,length(Ki))

    # run tested function
    Sj,Sj_bg,_,_ = simulateSpectrum(Fνj,hνj,Δνj,Ki,ΔKi,T,Kdisc,Be0,σ_cs_fg,σ_bg_vec)

    # results 
    (length(Sj)==length(Ki)) & (!isnan(sum(Sj))) & (!isinf(sum(Sj))) & (length(Sj_bg)==length(Ki)) & (!isnan(sum(Sj_bg))) & (!isinf(sum(Sj_bg)))
end



#######################################
##           run tests               ##
#######################################

@testset "Testing ATTIRE kernel" begin
    @test test_φi()
    @test test_Φi()
end
  
@testset "Testing ATTIRE light" begin
    @test test_sourceSpread()
end

@testset "Testing ATTIRE noise" begin
    @test test_countElectrons()
end

@testset "Testing ATTIRE analyzer" begin
    @test test_T_r4000()
    @test test_simulateSpectrum()
end