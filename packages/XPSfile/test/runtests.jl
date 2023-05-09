# XPSfile test bench

using Test
using XPSfile
using DataFrames
using Query




function test_dataAndFit_xlsx2df()
    # regular expression to navigate the files
    regData    = r"Eph ?= ?"                  # pattern found in the data sheet's name that does not appear in other sheets
    regFit     = r"[Ff]itt?ing(_| )?results?" # pattern found in fit sheet's name
    regEph     = r"=[0-9]*eV"                 # pattern for reading the photon energy
    sortBy     = Symbol("Photon energy")      # column name for the photon energy  (sort the data by increasing order of photon energy)
    thenSortBy = Symbol("Binding energy")     # colmun name for the binding energy (for the data sharing the same photon energy, sort them by binding energy)

    # load data from data_1
    dictAllData,df_fit,Ndata = dataAndFit_xlsx2df("data_1.xlsx"; regData=regData, regFit=regFit, regEph=regEph, sortBy=sortBy, thenSortBy=thenSortBy);

    # results 
    (length(dictAllData)==Ndata) & (typeof(df_fit)==DataFrames.DataFrame)
end

function test_dataAndFit_xlsx2df_missing()
    # regular expression to navigate the files
    regData    = r"Eph ?= ?"                  # pattern found in the data sheet's name that does not appear in other sheets
    regFit     = r"[Ff]itt?ing(_| )?results?" # pattern found in fit sheet's name
    regEph     = r"=[0-9]*eV"                 # pattern for reading the photon energy
    sortBy     = Symbol("Photon energy")      # column name for the photon energy  (sort the data by increasing order of photon energy)
    thenSortBy = Symbol("Binding energy")     # colmun name for the binding energy (for the data sharing the same photon energy, sort them by binding energy)

    # load data from data_2 (contains missing entries)
    dictAllData,df_fit,Ndata = dataAndFit_xlsx2df_missing("data_2.xlsx"; regData=regData, regFit=regFit, regEph=regEph, sortBy=sortBy, thenSortBy=thenSortBy);

    # results
    (length(dictAllData)==Ndata) & (typeof(df_fit)==DataFrames.DataFrame)
end





function test_dataAndMeta_xlsx2df()
    dictAllData,df,symbolDict = dataAndMeta_xlsx2df("data.xlsx",r"^hν_[0-9]*",r"meta");
    (length(dictAllData)==nrow(df)) & (length(symbolDict)==nrow(df))
end

function test_model_xlsx2df()
    dictAllGeom,symbolDictInv = model_xlsx2df("model.xlsx");    
    length(symbolDictInv)==length(dictAllGeom) 
end

@testset "XPSfile loading data and model" begin
    @test test_dataAndFit_xlsx2df()
    @test test_dataAndFit_xlsx2df_missing()
    @test test_dataAndMeta_xlsx2df()
    @test test_model_xlsx2df()
end




function test_curveFromFit()
    # Dataframe (fits)
    bind_sym   = Symbol("Binding energy");
    shift_sym  = Symbol("Peak shift");
    gauss_sym  = Symbol("FWHM(G)")
    loren_sym  = Symbol("FWHM(L)");
    area_sym   = Symbol("Area");
    dfPeaks = DataFrame(Dict(bind_sym=>[2.0; 4.0; 6.0], shift_sym=>[0.5; 0.2; 0.3], gauss_sym=>[0.6e3; 0.7e3; 0.5e3], loren_sym=>[0.1; 0.1; 0.1], area_sym=>[1.0; 2.0; 3.0]));
    Be       = collect(0.0:0.1:8.0);

    # curve from fits
    Npeaks, BePeak, σePeak, σePeakL, AePeak, σ_peak = curveFromFit(dfPeaks,Be,true;bind_sym=bind_sym, shift_sym=shift_sym, gauss_sym=gauss_sym, area_sym=area_sym, loren_sym=loren_sym);

    # results
    (Npeaks==3) & (BePeak==[2.0+0.5; 4.0+0.2; 6.0+0.3]) & (length(σePeak)==3) & (length(σePeakL)==3) & (AePeak==[1.0; 2.0; 3.0]) & (isapprox(sum(σ_peak)*(Be[2]-Be[1]),3.0;rtol=0.001))
end

@testset "XPSfile create curve from fit data" begin
    @test test_curveFromFit()
end