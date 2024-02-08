## 1-moment scheme

A scheme that predicts one moment of the particle size distribution (PSD).
Typically it's the 3rd moment of PSD.
The scheme's prognostic variable is the specific humidity (mass fraction)
  of water in each category, which is proportional to the 3rd moment of PSD.

## 2-moment scheme

A scheme that predicts two moments of the particle size distribution (PSD).
Typically those are the 3rd and the 0th moments of PSD.
The scheme's prognostic variables are the specific humidity (mass fraction)
  and number concentration of particles in each category,
  which are proportional to the 3rd and 0th moment of (PSD)

## Aerosol activation scheme

Aerosol particles serve as condensation nuclei for forming cloud droplets.
An aerosol activation scheme predicts the number concentrations
  of newly formed cloud droplets for a given population of aerosol particles.
The scheme is needed when using 2-moment microphysics.

## Ice nucleation scheme

Aerosol particles and cloud droplets serve as nuclei for forming ice crystals.
The pathways include water vapor deposition on dust, heterogeneous and
  homogeneous freezing of water droplets.
Ice nucleation schemes are needed to predict the number concentrations
  of newly formed ice crystals.
The schemes are needed when using 2-moment microphysics scheme.
