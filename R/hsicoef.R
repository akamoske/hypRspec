#' This function first reads in a corrected hdf5 file and applies a supplied coefficient to each band.
#' Coefficients may be from a PLSR, PCA, etc. Output is a raster file of same spatial resolution as input.
#'
#' @param hy.file hdf5 file containing hyperspectral imagery and associated metadata
#' @param metadata.path hdf5 path to reflectance metadata
#' @param reflectance.path hdf5 path to reflectance data
#' @param wavelength.path hdf5 path to wavelength metadata
#' @param coordinate.path hdf5 path to coordinate data
#' @param coef.csv csv containing coefficients
#' @param inter does the CSV contain an intercept and is it in the first position (TRUE/FALSE)
#' @param scale.data were the coefficients generated with scaled data and thus data needs to be scaled before applied coefs (TRUE/FALSE)
#' @param band.combo combination of bands index values that were used to generate the coefficients
#' @return A raster with all coefficients applied
#' @export
#' 
#' package requirement: rhdf5, terra (replacing "raster")

hsi.coef <- function(hy.file, metadata.path, reflectance.path, wavelength.path, 
                     coordinate.path, coef.csv, inter, scale.data, band.combo){
  
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
  
  # read in the coordinate infomation
  map.info <- h5read(file = hy.file,
                     name = coordinate.path)
  
  # save the crs projection data
  crs.proj <- base::paste0("+init=epsg:", map.info$`EPSG Code`)
  
  # pull out the map extent info
  map.info <- strsplit(map.info$Map_Info, split = ",", fixed = TRUE)
  x.min <- as.numeric(map.info[[1]][4])
  y.max <- as.numeric(map.info[[1]][5])
  
  if (inter == TRUE){
    
    # read in the PLSR coefficient data
    plsr.coef <- read.csv(coef.csv)
    
    # save the intercept 
    plsr.inter <- plsr.coef[1,2]
    
    # save just the coefficients
    plsr.coef <- plsr.coef[2:nrow(plsr.coef),]
    
    # set the matrix stack index
    s <- 1
    
    # first lets calculate the coefficients we will apply to all the images
    # if we are processing the entire image then we will remove the noisy bands
    # we can run this instead: c(25:194, 215:284, 325:403)
    # for just processing rbg images we can use: c(53,35,19)
    for (q in band.combo) {
      
      # lets read in the band and clean it up like we need before
      refl.array <- h5read(file = hy.file,
                           name = reflectance.path,
                           index = list(q, 1:n.cols, 1:n.rows))
      refl.matrix <- refl.array[1,,]
      refl.matrix[refl.matrix == data.ignore.val] <- NA
      refl.matrix <- refl.matrix / scale.fact.val
      
      # find the coef
      wl.coef <- plsr.coef[s,2]
      
      print(paste0("applying coefficient to ", plsr.coef[s,1]))
      
      # apply the coefficient
      if (scale.data == TRUE) {
        # plsr.matrix <- raster::scale(refl.matrix, center = TRUE, scale = FALSE) * wl.coef
        plsr.matrix <- scale(refl.matrix, center = TRUE, scale = FALSE) * wl.coef
      }
      if (scale.data == FALSE) {
        plsr.matrix <- refl.matrix * wl.coef
      }
      
      if (s == 1){
        # stack this raster
        hsi.matrix <- plsr.matrix
        hsi.matrix <- hsi.matrix + plsr.inter
      } else {
        hsi.matrix <- hsi.matrix + plsr.matrix
      }
      
      s <- s + 1
      
    } 
    
    #---------------------------------------------------------------------------------------------------
    # lets make a raster
    #---------------------------------------------------------------------------------------------------
    
    print(paste0("transforming flightline to a raster."))
    
    # convert the matrix to a raster
    # refl.raster <- raster(hsi.matrix, crs = crs.proj)
    refl.raster <- rast(hsi.matrix, crs = crs.proj)
    
    # we need to transpose the raster
    # refl.raster <- raster::t(refl.raster)
    refl.raster <- t(refl.raster)
    
    # find the dimensions of our raster
    y.dim <- dim(refl.raster)[1]
    x.dim <- dim(refl.raster)[2]
    
    # set the x.max and y.min
    x.max <- x.min + x.dim
    y.min <- y.max - y.dim
    
    # create an extent object
    # raster.ext <- extent(x.min, x.max, y.min, y.max)
    raster.ext <- ext(x.min, x.max, y.min, y.max)
    
    # set the spatial extent of the raster
    # extent(refl.raster) <- raster.ext
    ext(refl.raster) <- raster.ext
   
    return(refl.raster)
    
  } 
  
  if (inter == FALSE) {
    
    # read in the PLSR coefficient data
    plsr.coef <- read.csv(coef.csv)
    
    # set the matrix stack index
    s <- 1
    
    # first lets calculate the coefficients we will apply to all the images
    # if we are processing the entire image then we will remove the noisy bands
    # we can run this instead: c(25:194, 215:284, 325:403)
    # for just processing rbg images we can use: c(53,35,19)
    for (q in band.combo) {
      
      # lets read in the band and clean it up like we need before
      refl.array <- h5read(file = hy.file,
                           name = reflectance.path,
                           index = list(q, 1:n.cols, 1:n.rows))
      refl.matrix <- refl.array[1,,]
      refl.matrix[refl.matrix == data.ignore.val] <- NA
      refl.matrix <- refl.matrix / scale.fact.val
      
      # find the coef
      wl.coef <- plsr.coef[s,2]
      
      print(paste0("applying coefficient to ", plsr.coef[s,1]))
      
      # apply the coefficient
      if (scale.data == TRUE) {
        # plsr.matrix <- raster::scale(refl.matrix, center = TRUE, scale = FALSE) * wl.coef
        plsr.matrix <- scale(refl.matrix, center = TRUE, scale = FALSE) * wl.coef
      }
      if (scale.data == FALSE) {
        plsr.matrix <- refl.matrix * wl.coef
      }
      
      if (s == 1){
        # stack this raster
        hsi.matrix <- plsr.matrix
      } else {
        hsi.matrix <- hsi.matrix + plsr.matrix
      }
      
      s <- s + 1
      
    } 
    
    #---------------------------------------------------------------------------------------------------
    # lets make a raster
    #---------------------------------------------------------------------------------------------------
    
    print(paste0("transforming flightline to a raster."))
    
    # convert the matrix to a raster
    # refl.raster <- raster(hsi.matrix, crs = crs.proj)
    refl.raster <- rast(hsi.matrix, crs = crs.proj)
    
    # we need to transpose the raster
    # refl.raster <- raster::t(refl.raster)
    refl.raster <- t(refl.raster)
    
    # find the dimensions of our raster
    y.dim <- dim(refl.raster)[1]
    x.dim <- dim(refl.raster)[2]
    
    # set the x.max and y.min
    x.max <- x.min + x.dim
    y.min <- y.max - y.dim
    
    # create an extent object
    # raster.ext <- extent(x.min, x.max, y.min, y.max)
    raster.ext <- ext(x.min, x.max, y.min, y.max)
    
    # set the spatial extent of the raster
    # extent(refl.raster) <- raster.ext
    ext(refl.raster) <- raster.ext
    
    return(refl.raster)
  }
}
