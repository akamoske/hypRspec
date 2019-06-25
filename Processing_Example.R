# Install canopyLazR from GitHub
install_github("akamoske/hypRspec")

# Load the library
library(hypRspec)

# Calculate the NDVI mask
ndvi <- ndvi.mask(hy.file = "D:/BRDF_TESTING/TALL_HDF5/NEON_D08_TALL_DP1_20180429_190316_reflectance.h5",
                  metadata.path = "/TALL/Reflectance/Reflectance_Data",
                  reflectance.path = "/TALL/Reflectance/Reflectance_Data",
                  wavelength.path = "/TALL/Reflectance/Metadata/Spectral_Data/Wavelength",
                  red.nm = 674,
                  nir.nm = 830,
                  ndvi.threshold = 0.5)

# Calculate the brightness mask
brightness <- brightness.mask(hy.file = "D:/BRDF_TESTING/TALL_HDF5/NEON_D08_TALL_DP1_20180429_190316_reflectance.h5",
                              metadata.path = "/TALL/Reflectance/Reflectance_Data",
                              reflectance.path = "/TALL/Reflectance/Reflectance_Data",
                              wavelength.path = "/TALL/Reflectance/Metadata/Spectral_Data/Wavelength")

# Apply the topographic correction
hsi.raster <- hsi.correction(hy.file = "D:/BRDF_TESTING/TALL_HDF5/NEON_D08_TALL_DP1_20180429_190316_reflectance.h5",
                             ndvi.mask = ndvi,
                             brightness.mask = brightness,
                             band.combo = c(25:194, 215:284, 325:403),
                             metadata.path = "/TALL/Reflectance/Reflectance_Data",
                             reflectance.path = "/TALL/Reflectance/Reflectance_Data",
                             wavelength.path = "/TALL/Reflectance/Metadata/Spectral_Data/Wavelength",
                             solar.az.path = "/TALL/Reflectance/Metadata/Logs/Solar_Azimuth_Angle",
                             solar.zn.path = "/TALL/Reflectance/Metadata/Logs/Solar_Zenith_Angle",
                             slope.path = "/TALL/Reflectance/Metadata/Ancillary_Imagery/Slope",
                             aspect.path = "/TALL/Reflectance/Metadata/Ancillary_Imagery/Aspect",
                             sensor.az.path = "/TALL/Reflectance/Metadata/to-sensor_Azimuth_Angle",
                             sensor.zn.path = "/TALL/Reflectance/Metadata/to-sensor_Zenith_Angle",
                             coordinate.path = "/TALL/Reflectance/Metadata/Coordinate_System",
                             ross = "thick",
                             li = "dense",
                             raster.res = 10)



