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

# full_measurement_model.png

plot of some of the simulated measurements
