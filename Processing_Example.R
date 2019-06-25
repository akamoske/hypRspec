# Install canopyLazR from GitHub
install_github("akamoske/hypRspec")

# Load the library
library(hypRspec)

ndvi <- ndvi.mask(hy.file = "D:/BRDF_TESTING/TALL_HDF5/NEON_D08_TALL_DP1_20180429_190316_reflectance.h5",
                  metadata.path = "/TALL/Reflectance/Reflectance_Data",
                  reflectance.path = "/TALL/Reflectance/Reflectance_Data",
                  wavelength.path = "/TALL/Reflectance/Metadata/Spectral_Data/Wavelength",
                  red.nm = 674,
                  nir.nm = 830,
                  ndvi.threshold = 0.5)
