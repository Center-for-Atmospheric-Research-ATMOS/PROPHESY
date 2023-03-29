# XPS data simulation and inversion: PROPHESY
  this suite of packages should contain tools for the estimation of the concentration profiles of a  chemical species across a microjet (i.e. a very small stream) probed with X-rays. The data are the spectra of the emitted photoelectrons from the core subshell in the ground state (i.e. 1s orbital).

The suite is divided into four package for specific tasks:

- [XPSpack](packages/XPSpack/):         model of the photoelectron signal for several sample geometry 
- [XPSinv](packages/XPSinv/):           inversion methods from convex optimzation framework [6]
- [XPSsampling](packages/XPSsampling/): sampling methods for estimating the covariance of the posterior distribution and other distributions
- [XPSfile](packages/XPSfile/):         tools for loading data files



# Dependence

In the Project.toml files you'll find the dependences list hereafter

## XPSpack

- LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
- Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
- Interpolations = "a98d9a8b-a2ab-59e6-89dd-64a1c18fca59"
- Printf = "de0858da-6303-5e67-8744-51eddeeeb8d7"
- NMOpt = "ad0eafcd-0556-4673-9a02-4d69c7f573d5"

The package is an unregistered package that can be intalled from the github repository. In Julia RPEL, type:

`] add https://github.com/matthewozon/NMOpt`

## XPSinv

- LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
- Statistics = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"

## XPSfile

- XLSX = "fdbf4ff8-1666-58a4-91e7-1b58723a45e0"
- DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
- Query = "1a8c2f83-1ff3-5112-b086-8aa67b057ba1"
- DataValues = "e7dc6d0d-1eca-5fa6-8ad6-5aecde8b7ea5"

# Installation

## Julia

Julia can be downloaded from <https://julialang.org/downloads/>

## packages

For each package the instalation method is in their respective installation section ([XPSpack](packages/XPSpack/README.md), [XPSinv](packages/XPSinv/README.md), [XPSsampling](packages/XPSsampling/README.md) and [XPSfile](packages/XPSfile/README.md))

## tips

Once Julia is installed, if you do not install the package through the package manager (Pkg or ]), you need to tell Julia where to find the packages by updating the env variable `LOAD_PATH`. In the directory where Julia is installed, you can create the directory config (if it does not already exist) and add the file startup.jl (or edit it if it already exists). For instance, adding the `XPSpack` package located at `/path/to/package/XPSpack/`, you can add in the startup file the line:

`push!(LOAD_PATH,"/path/to/package/XPSpack/")`

where `/path/to/package/` needs to be change to your configuration.



On linux systems, the startup file is likely found at `~/.julia/config/startup.jl`.




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
- [12] M. Ozon, K. Tumashevich and N. L. Prisle , Quantitative alignment parameter estimation for analyzing X-ray photoelectron spectra, (submitted) Journal of Synchrotron Radiation, 2023, Vol. -, p. - ([DOI: 10.1107/-](https://www.doi.org/10.1107/-))
- [13]: M.Ozon, K. Tumashevich, J. J. Lin and N. L. Prisle , Inversion model for extracting chemically resolved depth profiles across liquid interfaces of various configurations from XPS data: PROPHESY, (submitted) Journal of Synchrotron Radiation, 2023, Vol. -, p. - ([DOI: 10.1107/-](https://www.doi.org/10.1107/-))
