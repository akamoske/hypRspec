#' THIS FUNCTION WILL OVERWRITE THE DATA THAT IS USED AS AN INPUT!!!!!!!!!!!
#' 
#' TO AVOID LOSING YOUR ORIGINAL DATA, MAKE SURE TO COPY THE ORIGINAL FILES TO ANOTHER LOCATION
#' BEFORE RUNNING THIS FUNCTION!!!!!!!!!
#' 
#' WE ARE NOT RESPONSIBLE FOR ANY LOST DATA FROM THE USE OF THIS FUNCTION.
#'
#' ---------------------------------------------------------------------------------------------------
#' 
#' Topographically and brdf correct a hdf5 file containing hyperspectral imagery and then write this data to a new hdf5 file.
#'
#' This function first reads in a hdf5 file and applies a topographic correction based on the following paper:
#' 
#' Maignan, F., F-M. Bréon, and R. Lacaze. "Bidirectional reflectance of Earth targets: Evaluation of analytical 
#' models using a large set of spaceborne measurements with emphasis on the Hot Spot." 
#' Remote Sensing of Environment 90.2 (2004): 210-220.(added)
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
#' Queally, Natalie, et al. "FlexBRDF: A flexible BRDF correction for grouped processing of airborne imaging spectroscopy flightlines." 
#' Journal of Geophysical Research: Biogeosciences 127.1 (2022): e2021JG006622. (added)
#' 
#' packages requirement: rhdf5
#'
#' @param hy.file hdf5 file containing hyperspectral imagery and associated metadata
#' @param ndvi.mask NDVI mask created with the ndvi.mask function
#' @param brightness.mask brightness mask created with the brightness.mask function
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
#' @param solar.zn.ref reference solar zenith angle (i.e., averaged SZA across flight linee; solar noon across a growing season, etc.)
#' @param sensor.zn.ref reference sensor zenith angle (nadir-viewing if VZA = 0)
#' @param relative.az.ref reference relative azimuth angle (doesn't matter if VZA = 0, it will be the same)
#' @param h_b_ratio h/b describes vegetation height (H: height-to-center-crown, B: crown vertical radius; R: crown horizontal radius)
#' @param b_r_ratop b/r describes crown shape
#' @return A matrix topographic and brdf corrected reflectance data (nadir-viewing, reference solar zenith angle)
#' @export
#' 

hsi.correct.write.2 <- function(hy.file, ndvi.mask, brightness.mask, 
                              metadata.path, reflectance.path, wavelength.path, 
                              solar.az.path, solar.zn.path, slope.path, aspect.path,
                              sensor.az.path, sensor.zn.path, coordinate.path, ross, 
                              li, solar.zn.ref, sensor.zn.ref, relative.az.ref, h_b_ratio, b_r_ratio, sample.bands){
  
  # get TileID
  file_string <- strsplit(hy.file, "_", fixed = "T")
  tileID <- paste0(file_string[[1]][5], "_", file_string[[1]][6])
  
  # lets look at the reflectance metadata
  refl.info <- h5readAttributes(hy.file, metadata.path)
  
  # lets read in the wavelength info
  wavelengths <- h5read(file = hy.file, 
                        name = wavelength.path)
  
  # lets save the dimensions of the dataset for future use
  n.rows <- refl.info$Dimensions[1]
  n.cols <- refl.info$Dimensions[2]
  n.bands <- refl.info$Dimensions[3]
  
  n.sample.bands <- length(sample.bands)
  sample.wavelengths <- wavelengths[sample.bands]
  
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
  solar.zn.ref <- (solar.zn.ref * pi)/180
  sensor.zn.ref <- (sensor.zn.ref *pi)/180
  relative.az.ref <- (relative.az.ref*pi)/180
  
  #---------------------------------------------------------------------------------------------------
  # lets calculate the topographic correction coefficients
  #---------------------------------------------------------------------------------------------------
  
  print("calculating topographic correction variables.")
  
  # Generate the cosine i
  # the cosine of the incidence angle (i), defined as the angle between the normal to the pixel 
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
  
  # first we need to use Eq 2. from Schlapfer et al. IEEE-TGARS 2015 which uses the inverse cosine (Maignan, 2004)
  phase <- acos(cos(solar.zn) * cos(sensor.zn) + sin(solar.zn) * sin(sensor.zn) * cos(relative.az))
  
  # for the Thick Ross Kernel - Eq 7. Wanner et al. JGRA 1995 (The original equation is "+" instead of "*"?) 
  # double check and revised based on https://www.star.nesdis.noaa.gov/pub/smcd/emb/promanov/OTHER/k_summ.pdf (M.c., 2024)
  ross.thick <- ((pi/2 - phase) * cos(phase) + sin(phase)) / (cos(sensor.zn) + cos(solar.zn)) - pi/4
  
  # for the Thin Ross Kernal - Eq 13. Wanner et al. JGRA 1995
  ross.thin <- ((pi/2 - phase) * cos(phase) + sin(phase)) / (cos(sensor.zn) * cos(solar.zn)) - pi/2
  
  #---------------------------------------------------------------------------------------------------
  # now we need to calculate the Ross Volumetric Scattering Kernel at nadir (sensor.zn.ref = 0)
  #---------------------------------------------------------------------------------------------------
  
  # first we need to use Eq 2. from Schlapfer et al. IEEE-TGARS 2015 which uses the inverse cosine
  phase.nad.sr <- acos(cos(solar.zn.ref) * cos(sensor.zn.ref) + sin(solar.zn.ref) * sin(sensor.zn.ref) * cos(relative.az.ref))
  
  # for the Thick Ross Kernel - Eq 13. Wanner et al. JGRA 1995 (Eq. 7? )
  ross.thick.nad <- ((pi/2 - phase.nad.sr) * cos(phase.nad.sr) + sin(phase.nad.sr)) / (cos(sensor.zn.ref) + cos(solar.zn.ref)) - pi/4
  
  # for the Thin Ross Kernel - Eq 13. Wanner et al. JGRA 1995
  ross.thin.nad <- ((pi/2 - phase.nad.sr) * cos(phase.nad.sr) + sin(phase.nad.sr)) / (cos(sensor.zn.ref) * cos(solar.zn.ref)) - pi/2
  
  #---------------------------------------------------------------------------------------------------
  # now we need to calculate the Li Geometric Scattering Kernel
  
  # we will use the constants from Colgan et al. 2012 in Remote Sensing (cited in the FlexBRDF paper as well, but they are for savanna trees)
  # h/b = 2 (vegetation height) ; b/r = 10 (crown shape)
  #---------------------------------------------------------------------------------------------------
  
  # first we need to implement Eq. 37,52. Wanner et al. JGRA 1995
  
  solar.zn.at <- atan(b_r_ratio * tan(solar.zn))
  sensor.zn.at <- atan(b_r_ratio * tan(sensor.zn))
  
  # next we need to use Eq 50. Wanner et al. JGRA 1995
  d <- sqrt((tan(solar.zn.at) ** 2) + (tan(sensor.zn.at) ** 2) - (2 * tan(solar.zn.at) * tan(sensor.zn.at) * cos(relative.az)))
  
  # next we need to use Eq 49. Wanner et al. JGRA 1995
  # we will need to restraint these values between -1 and 1 to not return NaN and Pi values
  # for more info see this stack exchange thread 
  # https://stackoverflow.com/questions/14026297/acos1-returns-nan-for-some-values-not-others
  
  t.num <- h_b_ratio * sqrt(d**2 + (tan(solar.zn.at) * tan(sensor.zn.at) * sin(relative.az)) ** 2)
  t.denom <- (1 / cos(solar.zn.at)) + (1 / cos(sensor.zn.at))
  t <- acos(pmax(pmin(h_b_ratio*(t.num/t.denom), 1.0), -1.0))
  
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
  solar.zn.at <- atan(b_r_ratio * tan(solar.zn.ref))
  sensor.zn.at <- atan(b_r_ratio * tan(sensor.zn.ref))
  
  # next we need to use Eq 50. Wanner et al. JGRA 1995
  d <- sqrt((tan(solar.zn.at) ** 2) + (tan(sensor.zn.at) ** 2) - (2 * tan(solar.zn.at) * tan(sensor.zn.at) * cos(relative.az.ref)))
  
  # next we need to use Eq 49. Wanner et al. JGRA 1995
  # we will need to restraint these values between -1 and 1 to not return NaN and Pi values
  # for more info see this stack exchange thread 
  # https://stackoverflow.com/questions/14026297/acos1-returns-nan-for-some-values-not-others
  
  t.num <- h_b_ratio * sqrt(d**2 + (tan(solar.zn.at) * tan(sensor.zn.at) * sin(relative.az.ref)) ** 2)
  t.denom <- (1 / cos(solar.zn.at)) + (1 / cos(sensor.zn.at))
  
  t <- acos(pmax(pmin((h_b_ratio*(t.num/t.denom)), 1.0), -1.0))
  
  # next we need to use Eq 33,48. Wanner et al. JGRA 1995
  o <- (1 / pi) * (t - sin(t) * cos(t)) * t.denom
  
  # next we need to use Eq 51. Wanner et al. JGRA 1995
  cos.phase <- cos(solar.zn.at) * cos(sensor.zn.at) + sin(solar.zn.at) * sin(sensor.zn.at) * cos(relative.az.ref)
  
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
  
  # memory clean up
  gc()
  rm(aspect)
  rm(cosine.i)
  rm(li)
  rm(li.dense)
  rm(li.dense.nad)
  rm(li.sparse)
  rm(li.sparse.nad)
  rm(rel.az)
  rm(ross)
  rm(ross.thick)
  rm(ross.thick.nad)
  rm(ross.thin)
  rm(ross.thin.nad)
  rm(slope)
  gc()
  
  coef.matrix <- matrix(data = NaN, nrow = n.sample.bands, ncol = 7)
  
  # lets correct the imagery and then overwrite it to the hdf5 file
  # for (q in 1:n.bands) {
  for (q in 1:n.sample.bands){
    
    bandID <- sample.bands[q]
    
    coef.matrix[q, 1] <- bandID
    coef.matrix[q, 2] <- wavelengths[q]
      
    print(paste0("applying topographic correction to band ", bandID, "."))
    
    # lets read in the band and clean it up like we need before
    refl.array <- h5read(file = hy.file,
                         name = reflectance.path,
                         index = list(bandID, 1:n.cols, 1:n.rows))
    refl.matrix <- refl.array[1,,]
    refl.matrix[refl.matrix == data.ignore.val] <- NA
    refl.matrix <- refl.matrix / scale.fact.val
    
    # png_filename <- paste0(hy.file, "_wavelength_", wavelengths[q], "nm.png")
    # png(png_filename)
    # image(refl.matrix)
    # dev.off()
    # 
    
    # memory clean up
    gc()
    rm(refl.array)
    gc()
    
    # lets apply the masks to this band
    refl.matrix <- ifelse(ndvi.mask, refl.matrix, NA)
    
    # lets transform the data into the appropriate shape for regression
    y <- matrix(refl.matrix, nrow = length(refl.matrix), ncol = 1)
    
    # lets run the regression now: # Eq 7. Soenen et al., IEEE TGARS 2005
    topo.lm <- lm(y ~ x.topo)
    print(topo.lm)
    
    # lets save the coefficients
    topo.coef <- topo.lm$coefficients
    coef.matrix[q, 3:4] <- topo.coef
    
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
    rm(refl.matrix)
    rm(cor.fact)
    gc()
    
    #---------------------------------------------------------------------------------------------------
    # lets apply the brdf correction to the topo corrected band
    #---------------------------------------------------------------------------------------------------
    
    print(paste0("applying brdf correction to band ", bandID, "."))
    
    # lets apply the brightness mask to this topo corrected band
    topo.matrix <- ifelse(brightness.mask, topo.cor, NA)
    
    # lets transform the data into the appropriate shape for regression
    y <- matrix(topo.matrix, nrow = length(topo.matrix), ncol = 1)
    
    # lets run the regression now
    brdf.lm <- lm(y ~ x.brdf)
    print(brdf.lm)
    
    # memory management
    gc()
    rm(y)
    gc()
    
    # lets save the coefficients
    brdf.coef <- brdf.lm$coefficients
    coef.matrix[q,5:7] <- brdf.coef
    
    # now lets apply the coefficients to the band - eq 5. Weyermann et al. IEEE-TGARS 2015
    brdf <- brdf.coef[[1]] + (li.mask * brdf.coef[[3]]) + (ross.mask * brdf.coef[[2]])
    brdf.nad <- brdf.coef[[1]] + (li.mask.n * brdf.coef[[3]]) + (ross.mask.n * brdf.coef[[2]])
    
    # lets find the correction factor: eq 4. Weyermann et al. IEEE-TGARS 2015
    brdf.cor <- brdf.nad / brdf
    
    # lets apply the correction factor to the band
    band.brdf <- topo.matrix * brdf.cor
    
    #---------------------------------------------------------------------------------------------------
    # lets clear up some memory
    #---------------------------------------------------------------------------------------------------
    
    gc()
    rm(topo.matrix)
    rm(brdf)
    rm(brdf.nad)
    rm(brdf.cor)
    rm(topo.cor)
    rm(brdf.lm)
    gc()
    
    #---------------------------------------------------------------------------------------------------
    # lets clean up the matrix and rewrite it to the hdf5 file
    #---------------------------------------------------------------------------------------------------
    
    # first we need to rescale the data so we don't have floating points
    band.brdf.scale <- round(scale.fact.val * band.brdf)
    
    # next we can replace the NA values with the data ignore value
    band.brdf.scale[is.na(band.brdf.scale)] <- data.ignore.val
    
    # next we can overwrite this reflectance data
    h5write(band.brdf.scale, 
            file = hy.file, 
            name = reflectance.path,
            index = list(bandID, NULL, NULL))
    
    # print an update
    print(paste0("Band #", bandID, " is finished!"))
    
  } 
  # save the coefficient matrix to table
  write.csv(coef.matrix, file = paste0(HSI_folder, "umbs_", tileID, "_correction_coef", ".csv"))
            
}
