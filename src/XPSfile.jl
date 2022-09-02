# push!(LOAD_PATH,"/home/mattoz/Dev/XPSinv.jl/src/")
# sudo nano /opt/julias/julia-1.3.1/etc/julia/startup.jl

"""

    little helper: loading data from files such as xlsx and formatting them into DataFrame and Dictionary

"""
module XPSfile

using XLSX
using DataFrames
using Query
using DataValues

export XLSX,DataFrames,Query,DataValues

# importing and formatting data from xlsx files to dtaframes/dictionaries
export dataAndFit_xlsx2df,curveFromFit

export dataAndMeta_xlsx2df, model_xlsx2df

export IGORcolumn_xlsx2df # quite specific 


include("XPSfile/XPSxlsx.jl")
end