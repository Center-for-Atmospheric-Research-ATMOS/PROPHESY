<!--[![XPSpack CI](https://github.com/Center-for-Atmospheric-Research-ATMOS/XPS-depth-inv/actions/workflows/CI_XPSpack.yml/badge.svg)](https://github.com/Center-for-Atmospheric-Research-ATMOS/XPS-depth-inv/actions/workflows/CI_XPSpack.yml)-->
<!--[![XPSpack CI](https://github.com/Center-for-Atmospheric-Research-ATMOS/XPS-depth-inv/actions/workflows/CI_XPSpack.yml/badge.svg?branch=beta)](https://github.com/Center-for-Atmospheric-Research-ATMOS/XPS-depth-inv/actions/workflows/CI_XPSpack.yml)-->

[![XPSpack CI](https://github.com/Center-for-Atmospheric-Research-ATMOS/XPS-depth-inv/actions/workflows/CI_XPSpack.yml/badge.svg)](https://github.com/Center-for-Atmospheric-Research-ATMOS/XPS-depth-inv/actions/workflows/CI_XPSpack.yml)
[![ATTIRE CI](https://github.com/Center-for-Atmospheric-Research-ATMOS/XPS-depth-inv/actions/workflows/CI_ATTIRE.yml/badge.svg)](https://github.com/Center-for-Atmospheric-Research-ATMOS/XPS-depth-inv/actions/workflows/CI_ATTIRE.yml)
[![XPSsampling CI](https://github.com/Center-for-Atmospheric-Research-ATMOS/XPS-depth-inv/actions/workflows/CI_XPSsampling.yml/badge.svg)](https://github.com/Center-for-Atmospheric-Research-ATMOS/XPS-depth-inv/actions/workflows/CI_XPSsampling.yml)
[![XPSinv CI](https://github.com/Center-for-Atmospheric-Research-ATMOS/XPS-depth-inv/actions/workflows/CI_XPSinv.yml/badge.svg)](https://github.com/Center-for-Atmospheric-Research-ATMOS/XPS-depth-inv/actions/workflows/CI_XPSinv.yml)
[![XPSfile CI](https://github.com/Center-for-Atmospheric-Research-ATMOS/XPS-depth-inv/actions/workflows/CI_XPSfile.yml/badge.svg)](https://github.com/Center-for-Atmospheric-Research-ATMOS/XPS-depth-inv/actions/workflows/CI_XPSfile.yml)



[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.8007333.svg)](https://doi.org/10.5281/zenodo.8007333)



# XPS data simulation and inversion: PROPHESY
  
  This suite of packages contains tools for the estimation of concentration profiles for chemical species across a microjet (i.e. a very small stream) probed with X-rays. The data are the spectra of the emitted photoelectrons from the core subshell in the ground state (i.e. 1s orbital).

The suite is divided into four package for specific tasks:

- [XPSpack](packages/XPSpack/):         model of the photoelectron signal for several sample geometry 
- [XPSinv](packages/XPSinv/):           inversion methods from convex optimzation framework [6]
- [XPSsampling](packages/XPSsampling/): sampling methods for estimating the covariance of the posterior distribution and other distributions
- [XPSfile](packages/XPSfile/):         tools for loading data files
- [ATTIRE](packages/ATTIRE/):           light model of the acquisition device (kinetic energy analyzer and light)


# Installation

## Julia

Julia can be downloaded from <https://julialang.org/downloads/>

## packages

For each package the instalation method is in their respective installation section ([XPSpack](packages/XPSpack/README.md), [XPSinv](packages/XPSinv/README.md), [XPSsampling](packages/XPSsampling/README.md), [XPSfile](packages/XPSfile/README.md), and [ATTIRE](packages/ATTIRE/README.md))

## tips

Once Julia is installed, if you do not install the package through the package manager (Pkg or ]), you need to tell Julia where to find the packages by updating the env variable `LOAD_PATH`. In the directory where Julia is installed, you can create the directory config (if it does not already exist) and add the file startup.jl (or edit it if it already exists). For instance, adding the `XPSpack` package located at `/path/to/package/XPSpack/`, you can add in the startup file the line:

`push!(LOAD_PATH,"/path/to/package/XPSpack/")`

where `/path/to/package/` needs to be change to your configuration.

On linux systems, the startup file is likely found at `~/.julia/config/startup.jl`.


# Example

## Spectra simulation 

![o1s_spectra](https://github.com/Center-for-Atmospheric-Research-ATMOS/PROPHESY/assets/7929598/922d2029-bf4a-4e09-a0fe-ba0e49110183)







## Alignment parameter estimation
![APE](https://github.com/Center-for-Atmospheric-Research-ATMOS/PROPHESY/assets/7929598/b3ccbf74-ff9d-4b15-96e9-ed8901ac76eb)

Alignment parameter simulation and estimation








## Sample emission model
![emission_model](https://github.com/Center-for-Atmospheric-Research-ATMOS/PROPHESY/assets/7929598/3a0a99d6-b83f-43e1-9665-7f71aab72368)

Geometry factor for a cylindrical sample [13]







## Depth profile reconstruction
![estimates_and_uncertainties_W10_model_err_smooth_edge_global_0003](https://github.com/Center-for-Atmospheric-Research-ATMOS/PROPHESY/assets/7929598/62651064-44a8-4b96-9dce-f2e85dc8396e)

Profile reconstruction in the case of $W_{10}$ for three levels of global attenuation length uncertainty [13]


# Refs

- [1] N. Ottosson et al., Photoelectron spectroscopy of liquid water and aqueous solution: Electron effective attenuation lengths and emission-angle anisotropy, Journal of Electron Spectroscopy and Related Phenomena, 2012, Vol. 177, p. 60 ([DOI: 10.1016/j.elspec.2009.08.007](https://www.doi.org/10.1016/j.elspec.2009.08.007))
- [2] D. Roy and D. Tremblay, Design of electron spectrometers, Reports on Progress in Physics, 1990, Vol. 53, No. 12, p. 1621 ([DOI: 10.1088/0034-4885/53/12/003](https://www.doi.org/10.1088/0034-4885/53/12/003))
- [3] S. Manson and J. Cooper, Photo-Ionization in the Soft x-Ray Range: 1 Z Dependence in a Central-Potential Model, Physical Review, 1968, Vol. 165, p. 126 ([DOI: 10.1103/PhysRev.165.126](https://www.doi.org/10.1103/PhysRev.165.126))
- [4] J. Yeh and I. Lindau, Atomic subshell photoionization cross sections and asymmetry parameters: 1⩽ Z⩽ 103, 1985, Vol. 32, No. 1, p. 1--155, ([DOI: 10.1016/0092-640X(85)90016-6](https://www.doi.org/10.1016/0092-640X\(85\)90016-6))
- [5] D. Stolzenburg, M. Ozon, M.  Kulmala, K. Lehtinen, K. Lehtipalo and J. Kangasluoma, 2022, Combining instrument inversions for sub-10 nm aerosol number size-distribution measurements, Journal of Aerosol Science, Vol. 159, p. 105862, [DOI: 10.1016/j.jaerosci.2021.105862](https://www.doi.org/10.1016/j.jaerosci.2021.105862)
- [6] A. Chambolle and T. Pock, 2011. A first-order primal-dual algorithm for convex problems with applications to imaging. Journal of mathematical imaging and vision, 40(1), pp.120-145, [DOI: 10.1007/s10851-010-0251-1](https://www.doi.org/10.1007/s10851-010-0251-1)
- [7] N. Metropolis et S. Ulam, The Monte Carlo method,  Journal of the American Statistical Association, vol. 44, no 247, 1949, p. 335–341, [DOI: 10.2307/2280232](https://www.doi.org/10.2307/2280232)
- [8] W.K. Hastings, Monte Carlo Sampling Methods Using Markov Chains and Their Applications, Biometrika, vol. 57, no 1, 1970, p. 97–109, [DOI: 10.2307/2334940](https://www.doi.org/10.2307/2334940)
- [9] S.-J. Baek, A. Park, Y.-J. Ahn and J. Choo,  Baseline correction using asymmetrically reweighted penalized least squares smoothing Analyst, Royal Society of Chemistry, 2015, Vol. 140, p. 250-257 [DOI: 10.1039/C4AN01061B](https://www.doi.org/10.1039/C4AN01061B)
- [10] A. P. Dempster, N. M. Laird  and D. B. Rubin,  Maximum likelihood from incomplete data via the EM algorithm, Journal of the royal statistical society: series B (methodological), Wiley Online Library, 1977, 39, 1-22 ([DOI: 10.1111/j.2517-6161.1977.tb01600.x](https://www.doi.org/10.1111/j.2517-6161.1977.tb01600.x))
- [11] S. Thürmer, R. Seidel, M. Faubel, W. Eberhardt, J. C. Hemminger, S. E. Bradforth and B. Winter, Photoelectron angular distributions from liquid water: Effects of electron scattering Physical review letters, Physical review letters, 2013, Vol. 111, p. 173005 ([DOI: 10.1103/PhysRevLett.111.173005](https://www.doi.org/10.1103/PhysRevLett.111.173005))
- [12] M. Ozon, K. Tumashevich and N. L. Prisle, Quantitative alignment parameter estimation for analyzing X-ray photoelectron spectra, Journal of Synchrotron Radiation, 2023, Vol. 30, p. 766-779 ([DOI: 
10.1107/S1600577523004150](https://www.doi.org/10.1107/S1600577523004150))
- [13] M. Ozon, K. Tumashevich, J. J. Lin and N. L. Prisle,  Inversion model for extracting chemically resolved depth profiles across liquid interfaces of various configurations from XPS data: PROPHESY, 2023, Vol. 30, p. - ([DOI: 10.1107/S1600577523006124](https://www.doi.org/10.1107/S1600577523006124))

