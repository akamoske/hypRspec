# hypRspec

R package to process airborne hyperspectral imagery.

## Information

For theory behind the functions in this package please see the citations in the individual scripts. 

For a package that performs similar functions in python, please see [HyTools.](https://github.com/EnSpec/HyTools-sandbox)  

# WARNING!!!!!!

The hsi.correct.write() function will overwrite the data that is used an input!

To avoid losing your original data, make sure to copy the original files to another location
before running this function.

The authors are not responsible for any lost data from the use of this function.

### Corresponding Author

Aaron G. Kamoske, PhD Candidate
   
  + [Michigan State University, Department of Geography, Environment, and Spatial Sciences](http://geo.msu.edu/)      
  + [ERSAM Lab](https://www.ersamlab.com/)   
  + akamoske@gmail.com

### Contributing Authors

Dr. Kyla M. Dahlin
  + [Michigan State University, Department of Geography, Environment, and Spatial Sciences](http://geo.msu.edu/)
  + [Michigan State University, Ecology, Evolutionary Biology, and Behavior Program](https://eebb.msu.edu/)
  + [ERSAM Lab](https://www.ersamlab.com/)
  + kdahlin@msu.edu

Meicheng Shen
  + [Michigan State University, Department of Geography, Environment, and Spatial Sciences](http://geo.msu.edu/)
  + [ERSAM Lab](https://www.ersamlab.com/)
  + shenmeic@msu.edu
  
## Installation

The easiest way to install `hypRspec` is via `install_github` from the `devtools` package:

```
# If you haven't already installed this package and its dependencies
install.packages("devtools")

# If you alread have devtools installed or just installed it
library(devtools)

# Install hypRspec from GitHub
install_github("akamoske/hypRspec")

# Load the library
library(hypRspec)

# If you need to install the rhdf5 package this is the best way to do it
source("http://bioconductor.org/biocLite.R")
biocLite("rhdf5")
```

## Example of usage (after installation)

Once the package is loaded into your R session, this is the an example of how to use the functions in this package
to process hyperspectral imagery:

```
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

# WARNING!!!
# THIS FUNCTION OVERWRITES THE INPUT FILE - PLEASE COPY THE ORIGINAL DATA BEFORE RUNNING THIS!!!

# Apply the corrections and OVERWRITE the hdf5 file
hsi.raster <- hsi.correct.write(hy.file = "D:/BRDF_TESTING/TALL_HDF5/NEON_D08_TALL_DP1_20180429_190316_reflectance.h5",
                                ndvi.mask = ndvi,
                                brightness.mask = brightness,
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
                                li = "dense")
                             
# extract reflectance data with a set of random points from the corrected HSI hdf5 file
hsi.refl.points <- hsi.random.extract("D:/BRDF_TESTING/TALL_HDF5/NEON_D08_TALL_DP1_20180429_190316_reflectance.h5", 
                                      metadata.path = "/TALL/Reflectance/Reflectance_Data",
                                      coordinate.path = "/TALL/Reflectance/Metadata/Coordinate_System",
                                      wavelength.path = "/TALL/Reflectance/Metadata/Spectral_Data/Wavelength",
                                      reflectance.path = "/TALL/Reflectance/Reflectance_Data",
                                      band.combo = c(25:194, 215:284, 325:403),
                                      number.pts = 370)
                                      
# extract reflectance data with a shapefile from the corrected HSI hdf5 file
hsi.refl.points <- hsi.shp.extract("D:/BRDF_TESTING/TALL_HDF5/NEON_D08_TALL_DP1_20180429_190316_reflectance.h5", 
                                   metadata.path = "/TALL/Reflectance/Reflectance_Data",
                                   coordinate.path = "/TALL/Reflectance/Metadata/Coordinate_System",
                                   wavelength.path = "/TALL/Reflectance/Metadata/Spectral_Data/Wavelength",
                                   reflectance.path = "/TALL/Reflectance/Reflectance_Data",
                                   band.combo = c(25:194, 215:284, 325:403),
                                   shp.file.loc = "./data/field_traits/shp/MLBS2018""
                                   shp.file.name = "MLBS2018_TOC_FoliarData_20190626")
                                      
# apply coefficients from a PCA or PLSR to the corrected HSI data and export a raster
hsi.coef.rst <- hsi.coef("D:/BRDF_TESTING/TALL_HDF5/NEON_D08_TALL_DP1_20180429_190316_reflectance.h5", 
                         metadata.path = "/TALL/Reflectance/Reflectance_Data",
                         coordinate.path = "/TALL/Reflectance/Metadata/Coordinate_System",
                         wavelength.path = "/TALL/Reflectance/Metadata/Spectral_Data/Wavelength",
                         reflectance.path = "/TALL/Reflectance/Reflectance_Data",
                         coef.csv = "./pc1_coef.csv",
                         inter = FALSE,
                         scale.data = TRUE,
                         band.combo = c(25:194, 215:284, 325:403))
                          
                             

```
## License

This project is licensed under the GNU GPUv2 License - see the [LICENSE.md](LICENSE.md) file for details

