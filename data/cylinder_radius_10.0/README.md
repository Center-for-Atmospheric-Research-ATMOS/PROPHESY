# data.xlsx

Each sheet in a data file corresponds to a photon energy.
The content of the data.xlsx file is labeled as follow:

  - Be: the binding energies
  - F:  photon flux
  - Ke: kinetic energies
  - Sbg: measured background signal (no noise)
  - Snoisy:      the total measured signal with noise
  - SpectrumA_1: the signal of interest without noise
  - Stot:        sum of the signal of interest and background without noise
  - T:           transmission gain
  - hν:          photon energy
  - λ:           penetration depth (attenuation depth, eal... or any other name)
  - μKe:         mean kinetic energy
  - σ_cs_dens:   cross section density is a probability density function that depends on the kinetic energy (integrate to 1, mind the dKe factor in the sum)
  - σ_tot:       total photo-ionization cross section values

For the simulation of water O1s spectra, the gas phase is also recorded:

  - SpectrumA_1_gas: the signal from the gas phase water
  - σ_cs_dens_gas:   cross section denity of the gas phase
  - α:               measured alignment parameter, i.e. the ratio of peak area of O1s liquid and gas phase

Meta-data appear in seperate sheets for the newest simulations. Some are already in the data, such as F, T, hν, λ, μKe and σ_tot, and other are new

  - peak_mode_gas and peak_mode_liq:   the binding energies for the gas pahse and liquid phase respectively
  - peak_width_gas and peak_width_liq: corss section peak widths for gas phase and liquid phase
  - peak_mode:                         binding energy of each peak
  - peak_probability:                  probability of each peak
  - peak_width:                        width of each peak



# model.xlsx

Each sheet in a data file corresponds to a prenetration depth.
The content of the model file is labeled

  - H:	       geometry factor assuming a uniform beam profile
  - H_true:	   geometry factor including the beam profile
  - hν:        photon energy
  - max_depth: maximum depth reach by the geomtry factor
  - model:     name of the model
  - r:         discretization radial distance (distance to the center of the sample in microns)
  - radius:    radius of the sharp edge cylinder (in microns)
  - x0:        coordinate of the analyzer aperture in microns (origin at the center of the sample)
  - y0:        coordinate of the analyzer aperture in microns (origin at the center of the sample)
  - z0:        coordinate of the analyzer aperture in microns (origin at the center of the sample)
  - δr:        distance outside the sharp edge volume (limit of the geometry factor at radius+δr)
  - λ:         penetration depth (attenuation depth, eal... or any other name)
  - ρ:         concentration (or density) profile

The model files potentially inculdes other parameters:

  - α:               true alignment parameter (true as in from the definition it was given in this framework)
  - H_gas:           approximated geometry factor (assuming uniform illumination) for the gas phase
  - H_gas_true:      geometry factor for the gas pahse (accounting for the beam profile)
  - r_gas:           discretization distance
  - xc, yc:          beam profile center (x axis is transversal and y axis is along the liquid microjet)
  - σx, σy:          beam profile spread
  - ρ_gas:           gas concentration profile (inversly proportional to the distance from the sample and depends on the pressure and temperature)

# full_measurement_model.png

plot of some of the simulated measurements
