# XPSpack

The package XPSpack is part of the PROPHESY suite. It is a package which covers some modelization aspects related to XPS experiments

 - geometry of the sample (plane, cylinder or sphere)
 - photon beam profile (e.g. Gaussian profile)
 - elemental total photoionization cross-section of C1s, O1s and S2p [1]
 - background removal [2]
 - photionization cross-section density estimation (model freee estimation)
 - peak fitting based on Expectation-Maximization algorithm [3]
 - attenuation length (form experimental data or reported in [4])
 - alignment parameter estimation [5]

This package was used and introduced in two manuscripts submitted to JSR [5,6]

# Dependence

XPSpack depends on some of the most common registered Julia packages, and one unregistered package, (see [Project.toml](Project.toml))
 - Statistics
 - LinearAlgebra
 - Interpolations
 - Printf
 - NMOpt is a package that implements some (Quasi) Newton methods (see [NMOpt](https://github.com/matthewozon/NMOpt))
 

# Installation

## Package manager

In the Julia REPL, type:

```
] add https://github.com/Center-for-Atmospheric-Research-ATMOS/XPS-depth-inv/tree/beta/packages/XPSpack
```

the package should be installed once Julia Pkg is done, however, if the installation fails, you may want to check the error/warning messages. If after fixing the bugs, the package still cannot be installed via Pkg, you may refer to the Manual section below.

## Manual

You can download the package and place it with you other unregistered packages, e.g. in your Dev folder.
You need to tell Julia where to find the packages by updating the env variable `LOAD_PATH`. In the directory where Julia is installed, you can create the directory config (if it does not already exist) and add the file startup.jl (or edit it if it already exists). For instance, adding the `XPSpack` package located at `/path/to/package/XPSpack/`, you can add in the startup file the line:

`push!(LOAD_PATH,"/path/to/package/XPSpack/")`

where `/path/to/package/` needs to be changed to your configuration.


On linux systems, the startup file is likely found at `~/.julia/config/startup.jl` or in the julia etc folder.



# Examples

In the [test](../../test/) folder, several examples of use of XPSpack are implemented:

 - [geometry/geom.jl](../../test/geometry/geom.jl): for several sample geometry, the geometry factors are compared
 - [outsideModel/geom_vapor.jl](../../test/outsideModel/geom_vapor.jl): show an example of geometry factor that takes into account the outside of the sample where the vapor or low-density-liquid exists and produces significant photoelectronic signal.
 - [alignmentParameter](../../test/alignmentParameter/): several examples of alignment parameter estimate are showcased (requires [XPSfile](../XPSfile/))
 - [bg_removal_and_proba_density.jl](../../test/bg_removal_and_proba_density.jl): illusatration of the background removal and photionization cross-section density estimation


# Refs

- [1] J. Yeh and I. Lindau, Atomic subshell photoionization cross sections and asymmetry parameters: 1⩽ Z⩽ 103, Atomic data and nuclear data tables, 1985, Vol. 32, No. 1, p. 1-155, ([DOI: 10.1016/0092-640X(85)90016-6](https://www.doi.org/10.1016/0092-640X\(85\)90016-6))
- [2] S.-J. Baek, A. Park, Y.-J. Ahn and J. Choo,  Baseline correction using asymmetrically reweighted penalized least squares smoothing Analyst, Royal Society of Chemistry, 2015, Vol. 140, p. 250-257 [DOI: 10.1039/C4AN01061B](https://www.doi.org/10.1039/C4AN01061B)
- [3] A. P. Dempster, N. M. Laird  and D. B. Rubin,  Maximum likelihood from incomplete data via the EM algorithm, Journal of the royal statistical society: series B (methodological), Wiley Online Library, 1977, 39, 1-22 ([DOI: 10.1111/j.2517-6161.1977.tb01600.x](https://www.doi.org/10.1111/j.2517-6161.1977.tb01600.x))
- [4] S. Thürmer, R. Seidel, M. Faubel, W. Eberhardt, J. C. Hemminger, S. E. Bradforth and B. Winter, Photoelectron angular distributions from liquid water: Effects of electron scattering Physical review letters, Physical review letters, 2013, Vol. 111, p. 173005 ([DOI: 10.1103/PhysRevLett.111.173005](https://www.doi.org/10.1103/PhysRevLett.111.173005))
- [5] M. Ozon, K. Tumashevich and N. L. Prisle , Quantitative alignment parameter estimation for analyzing X-ray photoelectron spectra, Journal of Synchrotron Radiation, 2023, Vol. -, p. - ([DOI: 10.1107/-](https://www.doi.org/10.1107/-))
- [6]: M.Ozon, K. Tumashevich, J. J. Lin and N. L. Prisle , Inversion model for extracting chemically resolved depth profiles across liquid interfaces of various configurations from XPS data: PROPHESY, Journal of Synchrotron Radiation, 2023, Vol. -, p. - ([DOI: 10.1107/-](https://www.doi.org/10.1107/-))
