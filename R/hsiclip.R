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
#' @param band.combo index of bands that will be processed
#' @param metadata.path hdf5 path to reflectance metadata
#' @param reflectance.path hdf5 path to reflectance data
#' @param wavelength.path hdf5 path to wavelength metadata
#' @param shp.file shape file containing point data to extract refelctance data
#' @return A matrix topographic and brdf corrected reflectance data
#' @export

hsi.clip <- function(hy.file, band.combo, 
                     metadata.path, reflectance.path, wavelength.path,
                     coordinate.path, shp.file.loc, shp.file.name){
  
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
  # now we can apply the topographic correction to each band
  #---------------------------------------------------------------------------------------------------
  
  # lets prevent R from writing the numbers in scientific format
  options(scipen = 999)
  
  # lets read in the shapefile
  toc.refl <- readOGR(shp.file.loc,
                      shp.file.name)
  
  # read in the coordinate infomation
  map.info <- h5read(file = hy.file,
                     name = coordinate.path)
  
  # save the crs projection data
  crs.proj <- base::paste0("+init=epsg:", map.info$`EPSG Code`)
  
  # pull out the map extent info
  map.info <- strsplit(map.info$Map_Info, split = ",", fixed = TRUE)
  x.min <- as.numeric(map.info[[1]][4])
  y.max <- as.numeric(map.info[[1]][5])
  
  # make an empty raster
  hsi.stack <- raster::stack()
  
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

    #---------------------------------------------------------------------------------------------------
    # lets make a raster
    #---------------------------------------------------------------------------------------------------
    
    print(paste0("extracting data from band ", q, "."))
    
    # convert the matrix to a raster
    refl.raster <- raster(refl.matrix, crs = crs.proj)
    
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
    
    #stack this raster
    hsi.stack <- stack(hsi.stack, refl.raster)
   
  } 
  
  return(hsi.stack)
  
}