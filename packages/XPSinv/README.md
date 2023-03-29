# XPSinv

The package XPSinv is part of the PROPHESY suite. It is a package which implements some convex optimization algorithms

 - Primal-dual [1]
 - Null Space Optimization [2]

This package was used and introduced in a manuscript submitted to JSR [3]

# Dependence

XPSinv depends on some of the most common registered Julia packages (see [Project.toml](Project.toml))
 - Statistics
 - LinearAlgebra


# Installation

## Package manager

In the Julia REPL, type:

```
] add https://github.com/Center-for-Atmospheric-Research-ATMOS/XPS-depth-inv/tree/beta/packages/XPSinv
```

the package should be installed once Julia Pkg is done, however, if the installation fails, you may want to check the error/warning messages. If after fixing the bugs, the package still cannot be installed via Pkg, you may refer to the Manual section below.

## Manual

You can download the package and place it with you other unregistered packages, e.g. in your Dev folder.
You need to tell Julia where to find the packages by updating the env variable `LOAD_PATH`. In the directory where Julia is installed, you can create the directory config (if it does not already exist) and add the file startup.jl (or edit it if it already exists). For instance, adding the `XPSinv` package located at `/path/to/package/XPSinv/`, you can add in the startup file the line:

`push!(LOAD_PATH,"/path/to/package/XPSinv/")`

where `/path/to/package/` needs to be changed to your configuration.


On linux systems, the startup file is likely found at `~/.julia/config/startup.jl` or in the julia etc folder.



# Examples

In the [test](../../test/) folder, several examples of use of XPSpack are implemented:

 - [reconstruction/estimation.jl](../../test/reconstruction/estimation.jl): concentration profile reconstruction from simulated data from different acquisition setup, e.g. number of photon energy acquisition (requires [XPSsampling](../XPSsampling/))
 - [reconstruction/profileEstimation.jl](../../test/reconstruction/profileEstimation.jl): reconstruction for one acquisition setup (requires [XPSsampling](../XPSsampling/))
 - [reconstruction/testTruncation.jl](../../test/reconstruction/testTruncation.jl): data simulation and reconstruction (requires [XPSsampling](../XPSsampling/) and [XPSpack](../XPSpack/))


# Refs

- [1] A. Chambolle and T. Pock, A first-order primal-dual algorithm for convex problems with applications to imaging, Journal of Mathematical Imaging and Vision, Springer, 2011, Vol. 40, p. 120-145, ([DOI: 10.1007/s10851-010-0251-1](https://www.doi.org/10.1007/s10851-010-0251-1))
- [2] D. Stolzenburg, M. Ozon, M. Kulmala, K. E. Lehtinen, K. Lehtipalo, and J. Kangasluoma, Combining instrument inversions for sub-10 nm aerosol number size-distribution measurements, Journal of Aerosol Science, Elsevier, 2022, Vol. 159, p. 105862, ([DOI: 10.1016/j.jaerosci.2021.105862](https://www.doi.org/10.1016/j.jaerosci.2021.105862))
- [3]: M.Ozon, K. Tumashevich, J. J. Lin and N. L. Prisle , Inversion model for extracting chemically resolved depth profiles across liquid interfaces of various configurations from XPS data: PROPHESY, (submitted) Journal of Synchrotron Radiation, 2023, Vol. -, p. - ([DOI: 10.1107/-](https://www.doi.org/10.1107/-))
