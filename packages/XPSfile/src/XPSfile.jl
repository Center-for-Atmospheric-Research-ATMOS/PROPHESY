#
# XPSfile.jl --
#
# XPSfile.jl is a helper module to load data from XLSX files
#
#------------------------------------------------------------------------------
#
# This file is part of the XPSfile module which is licensed under CC-BY 4.0 license.
#
# Copyright (C) 2022,  Matthew Ozon, Konstantin Tumashevich and Nønne L. Prisle.
#
#------------------------------------------------------------------------------


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
export dataAndFit_xlsx2df, dataAndFit_xlsx2df_missing,curveFromFit

export dataAndMeta_xlsx2df, model_xlsx2df

export IGORcolumn_xlsx2df # quite specific 


include("XPSxlsx.jl")
end
