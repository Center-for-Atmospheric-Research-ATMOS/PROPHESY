# XPSfile

The package XPSfile is part of the PROPHESY suite. This package intends to make dealing with loading and manipulating data files easier. It is highly specialized for XPS data files such as those exported from IGOR. It relies on querries, regular expressions and XLS file format.

# Dependence

XPSfile depends on some of registered Julia packages (see [Project.toml](Project.toml))

 - XLSX
 - DataFrames
 - Query
 - DataValues


# Installation

## Package manager

In the Julia REPL, type:

```
] add https://github.com/Center-for-Atmospheric-Research-ATMOS/XPS-depth-inv/tree/beta/packages/XPSfile
```

the package should be installed once Julia Pkg is done, however, if the installation fails, you may want to check the error/warning messages. If after fixing the bugs, the package still cannot be installed via Pkg, you may refer to the Manual section below.

## Manual

You can download the package and place it with you other unregistered packages, e.g. in your Dev folder.
You need to tell Julia where to find the packages by updating the env variable `LOAD_PATH`. In the directory where Julia is installed, you can create the directory config (if it does not already exist) and add the file startup.jl (or edit it if it already exists). For instance, adding the `XPSfile` package located at `/path/to/package/XPSfile/`, you can add in the startup file the line:

`push!(LOAD_PATH,"/path/to/package/XPSfile/")`

where `/path/to/package/` needs to be changed to your configuration.


On linux systems, the startup file is likely found at `~/.julia/config/startup.jl` or in the julia etc folder.




# Examples

In the [test](../../test/) folder, several examples of use of XPSfile are implemented:

 - [estimationFromDataFile](../../test/estimationFromDataFile/): concentration profile reconstruction from data files (requires [XPSpack](../XPSpack/), [XPSinv](../XPSinv/) and [XPSsampling](../XPSsampling/))
 - [alignmentParameter](../../test/alignmentParameter/): alignment parameter estimation from data files (requires [XPSpack](../XPSpack/))