#' Convert a hdf5 file containing hyperspectral imagery into a topographically corrected reflectance matrix.
#'
#' This function reads in a hdf5 file and create a topographically corrected reflectance matrix.
#' 
#' For theory behind this function please see:
#' 
#' Soenen, S.A., Peddle, D.R., and Coburn, C.A., 2005. SCS+C: A Modified Sun-Canopy-Sensor Topographic Correction 
#' in Forested Terrain. IEEE Transactions on Geoscience and Remote Sensing, 43(9): 2148-2159.
#'
#' @param hy.file hdf5 file containing hyperspectral imagery and associated metadata
#' @param ndvi.mask NDVI mask created with the ndvi.mask function
#' @param metadata.path hdf5 path to reflectance metadata
#' @param reflectance.path hdf5 path to reflectance data
#' @param wavelength.path hdf5 path to wavelength metadata
#' @param solar.az.path hdf5 path to solar azimuth data
#' @param solar.zn.path hdf5 path to solar zenith data
#' @param slope.path hdf5 path to slope data
#' @param aspect.path hdf5 path to aspect data
#' @param band.combo index of bands that will be processed
#' @return A list of topographically corrected matrices
#' @export

topo.correction <- function(hy.file, ndvi.mask, metadata.path, reflectance.path, wavelength.path, solar.az.path,
                            solar.zn.path, slope.path, aspect.path, band.combo){
  
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
  
  # now we will save several pieces of data
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

  # now we need to convert these to radians
  slope <- (slope * pi) / 180
  aspect <- (aspect * pi) / 180
  solar.az <- (solar.az * pi) / 180
  solar.zn <- (solar.zn * pi) / 180

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
  
  # create an empty list to store the corrected matrices in
  topo.list <- list()
  
  # set up our list index
  q <- 1
  
  # first lets calculate the coefficients we will apply to all the images
  # if we are processing the entire image then we will remove the noisy bands
  # we can run this instead: c(25:194, 215:284, 325:403)
  # for just processing rbg images we can use: c(53,35,19)
  for (i in c(band.combo)) {
    
    # do some updating
    print(paste0("Processing band ", i, "!"))
    
    # lets read in the band and clean it up like we need before
    refl.array <- h5read(file = hy.file,
                         name = reflectance.path,
                         index = list(i, 1:n.cols, 1:n.rows))
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
    
    # save our topographic corrected matrix to the list
    topo.list[[i]] <- topo.cor
    
    # update our list index
    q <- q + 1
   
  }
  
  # return the list of corrected matrices
  return(topo.list)
  
}