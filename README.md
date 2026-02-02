# Coregistration Automation for Drone Multispectral Imagery

# Overview

This repository presents a fully reproducible R-based workflow for automatic co-registration of drone multispectral imagery without the use of manual Ground Control Points (GCPs).

The method uses a correlation-based pixel-shift approach to align multispectral bands captured by multi-lens sensors, demonstrated using MicaSense RedEdge drone data.

The workflow includes:

 - Automatic band alignment

 - Saving aligned multispectral raster stacks

 - Visual assessment using false-color composites

# Motivation

Multispectral drone cameras capture different spectral bands using separate lenses, which introduces spatial misalignment between bands. Traditional correction methods rely on manual GCP selection, which is:

- Time-consuming
- Subjective
- Not scalable for large datasets

This project explores whether automatic band-to-band alignment can be achieved using simple correlation-based methods in R.

# Data Description
Sensor

MicaSense RedEdge

Bands Used: 
 - File Name	Band	Wavelength
 - IMG_0455_1.tif	Blue	475 nm
 - IMG_0455_2.tif	Green	560 nm
 - IMG_0455_3.tif	Red	668 nm
 - IMG_0455_4.tif	NIR	840 nm

# Methodology Summary

 - Raw multispectral bands are loaded as rasters.

 - Bands are converted to matrices for pixel-level processing.

 - The Green band is selected as the spatial reference.

 - Other bands (Blue, Red, NIR) are shifted pixel-by-pixel.

For each possible shift:

Overlapping pixels are compared

Pearson correlation is computed

The shift with the maximum correlation is selected.

Aligned bands are saved as a new raster stack.

Visual assessment is performed using false-color composites.

Why Green as Reference?

The Green band is used as the reference because:

 - It typically has high signal-to-noise ratio

 - It captures strong structural features (vegetation, edges)

 - It lies centrally in the visible spectrum

 - It correlates well with both visible and NIR bands

Repository Structure
coregistration_automation_drone_imagery/
│
├── coregistration_automation_feasibility.Rmd       # Reproducible R Markdown report
├── MicaSense_4Band_Aligned.tif   # Output aligned raster stack
├── IMG_0455_1.tif                # Blue band
├── IMG_0455_2.tif                # Green band
├── IMG_0455_3.tif                # Red band
├── IMG_0455_4.tif                # NIR band
└── README.md

Software Requirements

R 

R packages: terra


# Applicability to Other Drone Data

The method is not sensor-specific and can be applied to:

Other MicaSense datasets

Similar multi-lens multispectral cameras

Any dataset where bands share the same spatial resolution

Adjustments may be required for:

Different spatial resolutions

Larger misalignments

Radiometric differences


# Limitations

Assumes translational (shift-only) misalignment

Does not correct rotation or scale differences

Computational cost increases with larger search windows

Author

Iris Nana Obeng
(MSc. Global Change Ecology)

