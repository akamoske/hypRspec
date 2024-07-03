#' Extract refelctance data from a set of random points from a topographically and brdf corrected hdf5 file
#' containing hyperspectral imagery.
#'
#' @param hy.file hdf5 file containing corrected hyperspectral imagery
#' @param metadata.path hdf5 path to the metadata
#' @param reflectance.path hdf5 path to the reflectance data
#' @param coordinate.path hdf5 path to the coordinate information
#' @param wavelength.path hdf5 path to the wavelength information
#' @param band.combo combination of wavelengths you want to extract data from
#' @param number.pts number of random points you want to extract from the file
#' @return A matrix topographic and brdf corrected reflectance data
#' @export
#' 
#' package requirement: rhdf5, raster, tools
#' replace "raster" package functions with "terra" package
#' install.packages('installr')
#' install.Rtools(check = TRUE)

hsi.random.extract <- function(hy.file, metadata.path, reflectance.path, 
                               coordinate.path, wavelength.path, band.combo, number.pts){
  
  # lets look at the reflectance metadata
  refl.info <- h5readAttributes(hy.file, metadata.path)
  
  # lets save the dimensions of the dataset for future use
  n.rows <- refl.info$Dimensions[1]
  n.cols <- refl.info$Dimensions[2]
  n.bands <- refl.info$Dimensions[3]
  
  # lets save the scale factor and the data ignore value
  scale.fact.val <- refl.info$Scale_Factor
  data.ignore.val <- refl.info$Data_Ignore_Value
  
  # lets read in the wavelength info
  wavelengths <- h5read(file = hy.file, 
                        name = wavelength.path)
  
  # read in the coordinate infomation
  map.info <- h5read(file = hy.file,
                     name = coordinate.path)
  
  # pull out the map extent info
  map.inform <- strsplit(map.info$Map_Info, split = ",", fixed = TRUE)
  x.min <- as.numeric(map.inform[[1]][4])
  y.max <- as.numeric(map.inform[[1]][5])
  
  # save the crs projection data
  crs.proj <- base::paste0("+init=epsg:", map.info$`EPSG Code`)
  
  # lets prevent R from writing the numbers in scientific format
  options(scipen = 999)
  
  # we need to make an empty matrix to store all the reflectance data in
  ext.mat <- matrix(ncol = length(band.combo) + 2, 
                    nrow = number.pts + 1)
  
  # set the wavelength names
  wl.names <- paste0("nm", wavelengths[band.combo])
  
  # set the column names
  ext.mat[1,1] <- "ID"
  ext.mat[1,2] <- "Fline"
  ext.mat[1, 3:ncol(ext.mat)] <- wl.names
  
  # set the ID variable
  ext.mat[2:nrow(ext.mat), 1] <- as.vector(1:number.pts)
  
  # pull out the file name
  clean.name <- tools::file_path_sans_ext(hy.file)
  clean.name <- strsplit(clean.name, "/")[[1]][4]
  ext.mat[2:nrow(ext.mat), 2] <- strsplit(clean.name, "_")[[1]][3]
  
  # set the index for the matrix column
  r <- 3 
  
  # start a new count for the loop
  w <- 1
  
  # if we are processing the entire image then we will remove the noisy bands
  # we can run this instead: c(25:194, 215:284, 325:403)
  # for just processing rbg images we can use: c(53,35,19)
  
  # loop through the wavelengths, convert them to a raster, and extract the data
  for (q in band.combo) {
    
    # lets read in the band and clean it up like we need before
    refl.array <- h5read(file = hy.file,
                         name = reflectance.path,
                         index = list(q, 1:n.cols, 1:n.rows))
    refl.matrix <- refl.array[1,,]
    refl.matrix[refl.matrix == data.ignore.val] <- NA
    refl.matrix <- refl.matrix / scale.fact.val
    
    # memory clean up
    gc()
    rm(refl.array)
    gc()
    
    #---------------------------------------------------------------------------------------------------
    # lets make a raster
    #---------------------------------------------------------------------------------------------------
    
    print(paste0("extracting data from band ", q, "."))
    
    # convert the matrix to a raster
    # refl.raster <- raster(refl.matrix, crs = crs.proj)
    refl.raster <- rast(refl.matrix, crs = crs.proj)
    
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
    ext(refl.raster) <- raster.ext
    
    # the first time we need to create the random spatial dataframe to extract data from in subsequent loops
    if (w == 1){
      
      # lets extract the reflectance data
      # ref.toc <- sampleRandom(refl.raster, size = number.pts, na.rm = TRUE, sp = TRUE)
      ref.toc <- spatSample(refl.raster, size = number.pts, na.rm = TRUE, as.points = TRUE)
      
      # now we want to save the extracted data
      ref.data <- ref.toc@data
      
      # lets add this into the right part of the matrix
      ext.mat[2:nrow(ext.mat), r] <- as.vector(ref.data$layer)
    }
    
    else {
      
      # extract the reflectance data with the above spatial points dataframe
      # ref.toc@data[,w] <- raster::extract(refl.raster, ref.toc, method = "simple")
      ref.toc@data[,w] <- extract(refl.raster, ref.toc, fun= NULL, method = "simple")
      
      # now we want to save the extracted data
      ref.data <- ref.toc@data[,w]
      
      # lets add this into the right part of the matrix
      ext.mat[2:nrow(ext.mat), r] <- as.vector(ref.data)
    }
    
    # set the matrix index
    r <- r + 1
    
    #---------------------------------------------------------------------------------------------------
    # lets clear up some memory before making rasters
    #---------------------------------------------------------------------------------------------------
    
    gc()
    rm(refl.raster)
    gc()
    
    # update the count
    w <- w + 1
    
  } 
  
  return(ext.mat)
  
}