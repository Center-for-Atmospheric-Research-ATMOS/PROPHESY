# XPSsampling

The package XPSsampling is part of the PROPHESY suite. It is a package which implements the tools for stochastic sampling of the inversion model

 - communication mechanisms
 - covariance 
 - Metropolis-Hasting sampling algorithm [1,2]

The inversion model, based on the maximization of an a posteriori probability density, carries intrinsic uncertainty. The model is not infinitely sharp, so, the reconstruction lies within an acceptable distance from the true concentration profile. Furthermore, the measurement model is parametric, meaning that it depends on the value of parameters such as the attenuation length. From sample drawn from the a posteriori model, it is possible to approximate the marginal mean and covariance for some parameters. In a manuscript submitted to JSR [3] we show an example of marginalization over the attenuation length parameter and the consequences on the profile reconstructions.

# Dependence

XPSsampling is standalone, no external package are used  (see [Project.toml](Project.toml))

# Installation

## Package manager

In the Julia REPL, type:

```
] add https://github.com/Center-for-Atmospheric-Research-ATMOS/XPS-depth-inv/tree/beta/packages/XPSsampling
```

the package should be installed once Julia Pkg is done, however, if the installation fails, you may want to check the error/warning messages. If after fixing the bugs, the package still cannot be installed via Pkg, you may refer to the Manual section below.

## Manual

You can download the package and place it with you other unregistered packages, e.g. in your Dev folder.
You need to tell Julia where to find the packages by updating the env variable `LOAD_PATH`. In the directory where Julia is installed, you can create the directory config (if it does not already exist) and add the file startup.jl (or edit it if it already exists). For instance, adding the `XPSsampling` package located at `/path/to/package/XPSsampling/`, you can add in the startup file the line:

`push!(LOAD_PATH,"/path/to/package/XPSsampling/")`

where `/path/to/package/` needs to be changed to your configuration.


On linux systems, the startup file is likely found at `~/.julia/config/startup.jl` or in the julia etc folder.



# Examples

In the [test](../../test/) folder, several examples of use of XPSsampling are implemented:

 - [distributionSampling/samplingPosteriorData.jl](../../test/distributionSampling/samplingPosteriorData.jl): example of samples drawn from the posterior distribution and use of the sample to compute quantites such as covariance (requires [XPSinv](../XPSinv/) and [XPSpack](../XPSpack/))
 - [reconstruction/testTruncation.jl](../../test/reconstruction/testTruncation.jl): data simulation and reconstruction (requires [XPSinv](../XPSinv/) and [XPSpack](../XPSpack/))


# Refs

- [1] N. Metropolis et S. Ulam, The Monte Carlo method,  Journal of the American Statistical Association, vol. 44, no 247, 1949, p. 335–341, [DOI: 10.2307/2280232](https://www.doi.org/10.2307/2280232)
- [2] W.K. Hastings, Monte Carlo Sampling Methods Using Markov Chains and Their Applications, Biometrika, vol. 57, no 1, 1970, p. 97–109, [DOI: 10.2307/2334940](https://www.doi.org/10.2307/2334940)
- [3]: M.Ozon, K. Tumashevich, J. J. Lin and N. L. Prisle , Inversion model for extracting chemically resolved depth profiles across liquid interfaces of various configurations from XPS data: PROPHESY, Journal of Synchrotron Radiation, 2023, Vol. -, p. - ([DOI: 10.1107/-](https://www.doi.org/10.1107/-))
