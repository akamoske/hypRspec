#' Extract refelctance data from a topographically and brdf corrected hdf5 file containing hyperspectral imagery.
#'
#' This function first reads in a hdf5 file and applies a topographic correction based on the following paper:
#' 
#' Soenen, S.A., Peddle, D.R., and Coburn, C.A., 2005. SCS+C: A Modified Sun-Canopy-Sensor Topographic Correction 
#' in Forested Terrain. IEEE Transactions on Geoscience and Remote Sensing, 43(9): 2148-2159.
#' 
#' Next a brdf correction is applied based on the following papers:
#' 
#' Colgan, M.S., Baldeck, C.A., Feret, J.B., and Asner, G.P., 2012. Mapping savanna tree species at ecosystem scales 
#' using support vector machine classification and BRDF correction on airborne hyperspectral and LiDAR data.
#' Remote Sensing, 4(11): 3462-3480.  
#' 
#' Collings, S., Caccetta, P., Campbell, N., and Wu, X., 2010. Techniques for BRDF correction of hyperspectral mosaics. 
#' IEEE Transactions on Geoscience and Remote Sensing, 48(10): 3733-3746.
#' 
#' Schlapfer, D., Richter, R., and Feingersh, T., 2015. Operational BRDF effects correction for wide-field-of-view 
#' optical scanners (BREFCOR). IEEE Transactions on Geoscience and Remote Sensing, 53(4): 1855-1864.
#' 
#' Wanner, W., Li, X., and Strahler, A.H., 1995. On the derivation of kernels for kernel-driven models of 
#' bidirectional reflectance. Journal of Geophysical Research: Atmospheres, 100(D10): 21077-21089.
#' 
#' Weyermann, J., Kneubuhler, M., Schlapfer, D., and Schaepman, M.E., 2015. Minimizing Reflectance Anisotropy 
#' Effects in Airborne Spectroscopy Data Using Ross-Li Model Inversion With Continuous Field Land Cover Stratification. 
#' IEEE Transactions on Geoscience and Remote Sensing, 53(11): 5814-5823.
#' 
#' Last, refelctance data is extracted from this topographic and brdf corrected imagery
#'
#' @param hy.file hdf5 file containing hyperspectral imagery and associated metadata
#' @param ndvi.mask NDVI mask created with the ndvi.mask function
#' @param brightness.mask brightness mask created with the brightness.mask function
#' @param band.combo index of bands that will be processed
#' @param metadata.path hdf5 path to reflectance metadata
#' @param reflectance.path hdf5 path to reflectance data
#' @param wavelength.path hdf5 path to wavelength metadata
#' @param solar.az.path hdf5 path to solar azimuth data
#' @param solar.zn.path hdf5 path to solar zenith data
#' @param slope.path hdf5 path to slope data
#' @param aspect.path hdf5 path to aspect data
#' @param sensor.az.path hdf5 path to sensor azimuth data
#' @param sensor.zn.path hdf5 path to sensor zenith data
#' @param coordinate.path hdf5 path to coordinate data
#' @param ross set to either "thick" or "thin" based on ross kernal needed
#' @param li set to either "dense" or "sparse" based on ross kernal needed
#' @param shp.file shape file containing point data to extract refelctance data
#' @return A matrix topographic and brdf corrected reflectance data
#' @export

hsi.extract <- function(hy.file, ndvi.mask, brightness.mask, band.combo, 
                        metadata.path, reflectance.path, wavelength.path, 
                        solar.az.path, solar.zn.path, slope.path, aspect.path,
                        sensor.az.path, sensor.zn.path, coordinate.path, ross, 
                        li, shp.file.loc, shp.file.name){
  
  # lets look at the reflectance metadata
  refl.info <- h5readAttributes(hy.file, metadata.path)
  
  # lets read in the wavelength info
  wavelengths <- h5read(file = hy.file, 
                        name = wavelength.path)
  
  # lets save the dimensions of the dataset for future use
  n.rows <- refl.info$Dimensions[1]
  n.cols <- refl.info$Dimensions[2]
  n.bands <- refl.info$Dimensions[3]
  
  # lets save the scale factor and the data ignore value
  scale.fact.val <- refl.info$Scale_Factor
  data.ignore.val <- refl.info$Data_Ignore_Value
  
  #---------------------------------------------------------------------------------------------------
  # lets save all the information that will not need to be repeated for each band
  #---------------------------------------------------------------------------------------------------
  
  # now we will save several pieces of data
  sensor.az <- h5read(file = hy.file,
                      name = sensor.az.path)
  sensor.zn <- h5read(file = hy.file,
                      name = sensor.zn.path)
  solar.az <- h5read(file = hy.file,
                     name = solar.az.path)
  solar.zn <- h5read(file = hy.file,
                     name = solar.zn.path)
  slope <- h5read(file = hy.file,
                  name = slope.path)
  aspect <- h5read(file = hy.file,
                   name = aspect.path)
  
  # we need to remove the no data values from all of these
  slope[slope == data.ignore.val] <- NA
  aspect[aspect == data.ignore.val] <- NA
  solar.az[solar.az == data.ignore.val] <- NA
  solar.zn[solar.zn == data.ignore.val] <- NA
  sensor.az[sensor.az == data.ignore.val] <- NA
  sensor.zn[sensor.zn == data.ignore.val] <- NA
  
  # now we need to convert these to radians
  slope <- (slope * pi) / 180
  aspect <- (aspect * pi) / 180
  solar.az <- (solar.az * pi) / 180
  solar.zn <- (solar.zn * pi) / 180
  sensor.az <- (sensor.az * pi) / 180
  sensor.zn <- (sensor.zn * pi) / 180
 
  #---------------------------------------------------------------------------------------------------
  # lets calculate the topographic correction coefficients
  #---------------------------------------------------------------------------------------------------
  
  print("calculating topographic correction variables.")
  
  # Generate the cosine i
  # the cosine of the incidence angle (i ), defined as the angle between the normal to the pixel 
  # surface and the solar zenith direction
  rel.az <- aspect - solar.az
  cosine.i <- cos(solar.zn) * cos(slope) + sin(solar.zn) * sin(slope) * cos(rel.az)
  
  # lets apply the NDVI mask to the cosine.i
  cosine.i.mask <- ifelse(ndvi.mask, cosine.i, NA)
  
  # lets reshape this for regression
  x.topo <- matrix(cosine.i.mask, nrow = length(cosine.i.mask), ncol = 1)
  
  # Eq 11. Soenen et al., IEEE TGARS 2005
  # cos(alpha)* cos(theta)
  # alpha -- slope (slope), theta -- solar zenith angle (solar_zn)
  c1 <- cos(solar.zn) * cos(slope)
  
  #---------------------------------------------------------------------------------------------------
  # now we need to calculate the Ross Volumetric Scattering Kernel
  #---------------------------------------------------------------------------------------------------
  
  print("calculating brdf correction variables.")
  
  # we need to calculate the relative azimuth first
  relative.az <- sensor.az - solar.az
  
  # first we need to use Eq 2. from Schlapfer et al. IEEE-TGARS 2015 which uses the inverse cosine
  phase <- acos(cos(solar.zn) * cos(sensor.zn) + sin(solar.zn) * sin(sensor.zn) * cos(relative.az))
  
  # for the Thick Ross Kernel - Eq 7. Wanner et al. JGRA 1995
  ross.thick <- ((pi/2 - phase) * cos(phase) + sin(phase)) / (cos(sensor.zn) * cos(solar.zn)) - pi/4
  
  # for the Thin Ross Kernal - Eq 13. Wanner et al. JGRA 1995
  ross.thin <- ((pi/2 - phase) * cos(phase) + sin(phase)) / (cos(sensor.zn) * cos(solar.zn)) - pi/2
  
  #---------------------------------------------------------------------------------------------------
  # now we need to calculate the Ross Volumetric Scattering Kernel at nadir (sensor.zn = 0)
  #---------------------------------------------------------------------------------------------------
  
  # first we need to use Eq 2. from Schlapfer et al. IEEE-TGARS 2015 which uses the inverse cosine
  phase.nad <- acos(cos(solar.zn) * cos(0) + sin(solar.zn) * sin(0) * cos(relative.az))
  
  # for the Thick Ross Kernel - Eq 13. Wanner et al. JGRA 1995
  ross.thick.nad <- ((pi/2 - phase.nad) * cos(phase.nad) + sin(phase.nad)) / (cos(0) * cos(solar.zn)) - pi/4
  
  # for the Thin Ross Kernal - Eq 13. Wanner et al. JGRA 1995
  ross.thin.nad <- ((pi/2 - phase.nad) * cos(phase.nad) + sin(phase.nad)) / (cos(0) * cos(solar.zn)) - pi/2
  
  #---------------------------------------------------------------------------------------------------
  # now we need to calculate the Li Geometric Scattering Kernel
  
  # we will use the constants from Colgan et al. 2012 in Remote Sensing
  # h/b = 2 ; b/r = 10
  #---------------------------------------------------------------------------------------------------
  
  # first we need to implement Eq. 37,52. Wanner et al. JGRA 1995
  solar.zn.at <- atan(10 * tan(solar.zn))
  sensor.zn.at <- atan(10 * tan(sensor.zn))
  
  # next we need to use Eq 50. Wanner et al. JGRA 1995
  d <- sqrt((tan(solar.zn.at) ** 2) + (tan(sensor.zn.at) ** 2) - (2 * tan(solar.zn.at) * tan(sensor.zn.at) * cos(relative.az)))
  
  # next we need to use Eq 49. Wanner et al. JGRA 1995
  # we will need to restraint these values between -1 and 1 to not return NaN and Pi values
  # for more info see this stack exchange thread 
  # https://stackoverflow.com/questions/14026297/acos1-returns-nan-for-some-values-not-others
  
  t.num <- 2 * sqrt(d**2 + (tan(solar.zn.at) * tan(sensor.zn.at) * sin(relative.az)) ** 2)
  t.denom <- (1 / cos(solar.zn.at)) + (1 / cos(sensor.zn.at))
  t <- acos(pmax(pmin((t.num/t.denom), 1.0), -1.0))
  
  # next we need to use Eq 33,48. Wanner et al. JGRA 1995
  o <- (1 / pi) * (t - sin(t) * cos(t)) * t.denom
  
  # next we need to use Eq 51. Wanner et al. JGRA 1995
  cos.phase <- cos(solar.zn.at) * cos(sensor.zn.at) + sin(solar.zn.at) * sin(sensor.zn.at) * cos(relative.az)
  
  # for the Sparse Li Kernel - Eq 32. Wanner et al. JGRA 1995
  li.sparse <- o - (1 / cos(solar.zn.at)) - (1 / cos(sensor.zn.at)) + 0.5 * (1 + cos.phase) * (1 / cos(sensor.zn.at))
  
  # for the Dense Li Kernel - Eq 47. Wanner et al. JGRA 1995
  li.dense <- (((1 + cos.phase) * (1 / cos(sensor.zn.at))) / (t.denom - o)) - 2
  
  #---------------------------------------------------------------------------------------------------
  # now we need to calculate the Li Geometric Scattering Kernel at nadir (sensor.zn = 0)
  
  # we will use the constants from Colgan et al. 2012 in Remote Sensing
  # h/b = 2 ; b/r = 10
  #---------------------------------------------------------------------------------------------------
  
  # first we need to implement Eq. 37,52. Wanner et al. JGRA 1995
  solar.zn.at <- atan(10 * tan(solar.zn))
  sensor.zn.at <- atan(10 * tan(0))
  
  # next we need to use Eq 50. Wanner et al. JGRA 1995
  d <- sqrt((tan(solar.zn.at) ** 2) + (tan(sensor.zn.at) ** 2) - (2 * tan(solar.zn.at) * tan(sensor.zn.at) * cos(relative.az)))
  
  # next we need to use Eq 49. Wanner et al. JGRA 1995
  # we will need to restraint these values between -1 and 1 to not return NaN and Pi values
  # for more info see this stack exchange thread 
  # https://stackoverflow.com/questions/14026297/acos1-returns-nan-for-some-values-not-others
  
  t.num <- sqrt(d**2 + (tan(solar.zn.at) * tan(sensor.zn.at) * sin(relative.az) ** 2))
  t.denom <- (1 / cos(solar.zn.at)) + (1 / cos(sensor.zn.at))
  t <- acos(pmax(pmin((2*(t.num/t.denom)), 1.0), -1.0))
  
  # next we need to use Eq 33,48. Wanner et al. JGRA 1995
  o <- (1 / pi) * (t - sin(t) * cos(t)) * t.denom
  
  # next we need to use Eq 51. Wanner et al. JGRA 1995
  cos.phase <- cos(solar.zn.at) * cos(sensor.zn.at) + sin(solar.zn.at) * sin(sensor.zn.at) * cos(relative.az)
  
  # for the Sparse Li Kernel - Eq 32. Wanner et al. JGRA 1995
  li.sparse.nad <- o - (1 / cos(solar.zn.at)) - (1 / cos(sensor.zn.at)) + 0.5 * (1 + cos.phase) * (1 / cos(sensor.zn.at))
  
  # for the Dense Li Kernel - Eq 47. Wanner et al. JGRA 1995
  li.dense.nad <- (((1 + cos.phase) * (1 / cos(sensor.zn.at))) / (t.denom - o)) - 2
  
  #---------------------------------------------------------------------------------------------------
  # lets clear up some memory since the only things we need from above are the kernels
  #---------------------------------------------------------------------------------------------------
  
  gc()
  rm(cos.phase)
  rm(d)
  rm(o)
  rm(phase)
  rm(phase.nad)
  rm(relative.az)
  rm(sensor.az)
  rm(sensor.zn)
  rm(t)
  rm(t.num)
  gc()
  
  #---------------------------------------------------------------------------------------------------
  # lets create our vegetation masks for the kernels we will use and then rearrange the data
  # for regression
  #---------------------------------------------------------------------------------------------------
  
  # lets apply the NDVI and brightness masks to the kernals based on which kernals we want to use
  if (ross == "thick"){
    print("Using thick ross.")
    ross.mask <- ifelse(ndvi.mask, ross.thick, NA)
    ross.mask <- ifelse(brightness.mask, ross.mask, NA)
    ross.mask.n <- ifelse(ndvi.mask, ross.thick.nad, NA)
    ross.mask.n <- ifelse(brightness.mask, ross.mask.n, NA)
  } 
  
  if (ross == "thin"){
    print("Using thin ross.")
    ross.mask <- ifelse(ndvi.mask, ross.thin, NA)
    ross.mask <- ifelse(brightness.mask, ross.mask, NA)
    ross.mask.n <- ifelse(ndvi.mask, ross.thin.nad, NA)
    ross.mask.n <- ifelse(brightness.mask, ross.mask.n, NA)
  }
  
  if(li == "dense"){
    print("Using dense li.")
    li.mask <- ifelse(ndvi.mask, li.dense, NA)
    li.mask <- ifelse(brightness.mask, li.mask, NA)
    li.mask.n <- ifelse(ndvi.mask, li.dense.nad, NA)
    li.mask.n <- ifelse(brightness.mask, li.mask.n, NA)
  }
  
  if(li == "sparse"){
    print("Using sparse li.")
    li.mask <- ifelse(ndvi.mask, li.sparse, NA)
    li.mask <- ifelse(brightness.mask, li.mask, NA)
    li.mask.n <- ifelse(ndvi.mask, li.sparse.nad, NA)
    li.mask.n <- ifelse(brightness.mask, li.mask.n, NA)
  }
  
  # lets transform the data into the appropriate shape for regression
  ross <- matrix(ross.mask, nrow = length(ross.mask), ncol = 1)
  li <- matrix(li.mask, nrow = length(li.mask), ncol = 1)
  
  # lets column bind the data we need together
  x.brdf <- cbind(ross, li)
  
  #---------------------------------------------------------------------------------------------------
  # now we can apply the topographic correction to each band
  #---------------------------------------------------------------------------------------------------
  
  # lets prevent R from writing the numbers in scientific format
  options(scipen = 999)
  
  # lets read in the shapefile
  toc.refl <- readOGR(shp.file.loc,
                      shp.file.name)
  
  # we need to make an empty matrix to store all the reflectance data in
  ext.mat <- matrix(ncol = length(band.combo) + 2, 
                    nrow = nrow(toc.refl) + 1)
  
  # set the wavelength names
  wl.names <- paste0("nm", wavelengths[band.combo])
  
  # set the column names
  ext.mat[1,1] <- "ID"
  ext.mat[1,2] <- "Fline"
  ext.mat[1, 3:ncol(ext.mat)] <- wl.names
  
  # make sure that the column variable is correctly named
  colnames(toc.refl@data) <- "ID"
  
  # set the ID variable
  ext.mat[2:nrow(ext.mat), 1] <- as.vector(toc.refl@data$ID)
  
  # pull out the file name
  ext.mat[2:nrow(ext.mat), 2] <- strsplit(hy.file, "_")[[1]][8]
  
  # set the index for the matrix column
  r <- 3 
  
  # read in the coordinate infomation
  map.info <- h5read(file = hy.file,
                     name = coordinate.path)
  
  # save the crs projection data
  crs.proj <- base::paste0("+init=epsg:", map.info$`EPSG Code`)
  
  # pull out the map extent info
  map.info <- strsplit(map.info$Map_Info, split = ",", fixed = TRUE)
  x.min <- as.numeric(map.info[[1]][4])
  y.max <- as.numeric(map.info[[1]][5])
  
  # first lets calculate the coefficients we will apply to all the images
  # if we are processing the entire image then we will remove the noisy bands
  # we can run this instead: c(25:194, 215:284, 325:403)
  # for just processing rbg images we can use: c(53,35,19)
  for (q in band.combo) {
    
    print(paste0("applying topographic correction to band ", q, "."))
    
    # lets read in the band and clean it up like we need before
    refl.array <- h5read(file = hy.file,
                         name = reflectance.path,
                         index = list(q, 1:n.cols, 1:n.rows))
    refl.matrix <- refl.array[1,,]
    refl.matrix[refl.matrix == data.ignore.val] <- NA
    refl.matrix <- refl.matrix / scale.fact.val
    
    # lets apply the masks to this band
    refl.matrix <- ifelse(ndvi.mask, refl.matrix, NA)
    
    # lets transform the data into the appropriate shape for regression
    y <- matrix(refl.matrix, nrow = length(refl.matrix), ncol = 1)
    
    # lets run the regression now: # Eq 7. Soenen et al., IEEE TGARS 2005
    topo.lm <- lm(y ~ x.topo)
    
    # lets save the coefficients
    topo.coef <- topo.lm$coefficients
    
    # Eq 8. Soenen et al., IEEE TGARS 2005
    if (topo.coef[[2]] == 0){
      c <- 0
    } else {
      c <- topo.coef[[1]] / topo.coef[[2]]
    }
    
    # find correction factor - Eq 11. Soenen et al., IEEE TGARS 2005
    cor.fact <- (c1 + c) / (cosine.i.mask + c)
    
    # apply the correction factor
    topo.cor <- refl.matrix * cor.fact
    
    #---------------------------------------------------------------------------------------------------
    # lets clear up some memory before applying the brdf correction
    #---------------------------------------------------------------------------------------------------
    
    gc()
    rm(topo.lm)
    gc()
    
    #---------------------------------------------------------------------------------------------------
    # lets apply the brdf correction to the topo corrected band
    #---------------------------------------------------------------------------------------------------
    
    print(paste0("applying brdf correction to band ", q, "."))
    
    # lets apply the brightness mask to this topo corrected band
    topo.matrix <- ifelse(brightness.mask, topo.cor, NA)
    
    # lets transform the data into the appropriate shape for regression
    y <- matrix(topo.matrix, nrow = length(topo.matrix), ncol = 1)
    
    # lets run the regression now
    brdf.lm <- lm(y ~ x.brdf)
    
    # lets save the coefficients
    brdf.coef <- brdf.lm$coefficients
    
    # now lets apply the coefficients to the band - eq 5. Weyermann et al. IEEE-TGARS 2015
    brdf <- brdf.coef[[1]] + (li.mask * brdf.coef[[3]]) + (ross.mask * brdf.coef[[2]])
    brdf.nad <- brdf.coef[[1]] + (li.mask.n * brdf.coef[[3]]) + (ross.mask.n * brdf.coef[[2]])
    
    # lets find the correction factor: eq 4. Weyermann et al. IEEE-TGARS 2015
    brdf.cor <- brdf.nad / brdf
    
    # lets apply the correction factor to the band
    band.brdf <- topo.matrix * brdf.cor
    
    #---------------------------------------------------------------------------------------------------
    # lets clear up some memory before applying the brdf correction
    #---------------------------------------------------------------------------------------------------
    
    gc()
    rm(brdf.lm)
    gc()
    
    #---------------------------------------------------------------------------------------------------
    # lets make a raster
    #---------------------------------------------------------------------------------------------------
    
    print(paste0("extracting data from band ", q, "."))
   
    # convert the matrix to a raster
    refl.raster <- raster(band.brdf, crs = crs.proj)
    
    # we need to transpose the raster
    refl.raster <- raster::t(refl.raster)
    
    # find the dimensions of our raster
    y.dim <- dim(refl.raster)[1]
    x.dim <- dim(refl.raster)[2]
    
    # set the x.max and y.min
    x.max <- x.min + x.dim
    y.min <- y.max - y.dim
    
    # create an extent object
    raster.ext <- extent(x.min, x.max, y.min, y.max)
    
    # set the spatial extent of the raster
    extent(refl.raster) <- raster.ext
    
    # lets extract the reflectance data
    ref.toc <- extract(x = refl.raster,
                       y = toc.refl,
                       method = "simple")
    
    # lets add this into the right part of the matrix
    ext.mat[2:nrow(ext.mat), r] <- ref.toc
    
    # set the matrix index
    r <- r + 1
    
    #---------------------------------------------------------------------------------------------------
    # lets clear up some memory before making rasters
    #---------------------------------------------------------------------------------------------------
    
    gc()
    rm(refl.raster)
    gc()
    
  } 
  
  return(ext.mat)
  
}