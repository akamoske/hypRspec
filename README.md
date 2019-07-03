# hypRspec

R package to process airborne hyperspectral imagery.

## Information

For theory behind the functions in this package please see the citations in the individual scripts. 

For a python version of these functions, please see [HyTools.](https://github.com/EnSpec/HyTools-sandbox)   

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

# Apply the corrections and export a raster
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
                             
# Apply the corrections and extract the reflectance data
hsi.refl <- hsi.extract(hy.file = "D:/Tests/BRDF_TESTING/TALL_HDF5/NEON_D08_TALL_DP1_20180429_190316_reflectance.h5",
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
                        shp.file.loc = "C:/PROCESSED_DATA/SHP_FILES/FOLIAR_DATA/TALL2018",
                        shp.file.name = "TALL2018_TOC_FoliarData_20190626")
                        
# Apply the PLSR coefficients and export a raster
hsi.plsr <- hsi.plsr(hy.file = hsi.files[i],
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
                     li = "dense",
                     plsr.csv = "./PROCESSED_DATA/PLSR_TOC/leaf.N_PLSR_Coefficients_6comp.csv"))
```
## License

This project is licensed under the GNU GPUv2 License - see the [LICENSE.md](LICENSE.md) file for details

