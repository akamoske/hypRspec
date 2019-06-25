#' Convert a hdf5 file containing hyperspectral imagery into a NDVI mask for use in later corrections.
#'
#' This function reads in a hdf5 file and create a NDVI matrix for use in the topographic and BRDF corrections that are
#' applied in later steps.
#' 
#' For theory behind this function please see:
#' 
#' Asner, G.P., Martin, R.E., Anderson, C.B., and Knapp, D.E., 2015. Quantifying forest canopy traits: Imaging spectroscopy
#' versus field survey. Remote Sensing of Environment, 158: 15-27.
#' 
#' Dahlin, K.M., Asner, G.P., Field., C.B., 2014. Linking vegetation patterns to environmental gradients and human impacts 
#' in a mediterranean-type island ecosystem. Landscape Ecology, 29(9): 1571-1585.
#'
#' @param hy.file hdf5 file containing hyperspectral imagery and associated metadata
#' @param metadata.path hdf5 path to reflectance metadata
#' @param reflectance.path hdf5 path to reflectance data
#' @param wavelength.path hdf5 path to wavelength metadata
#' @param red.nm wavelength (nm) to use for red band
#' @param nir.nm wavelength (nm) to use for nir band
#' @param ndvi.threshold threshold for NDVI mask
#' @return A matrix containing the NDVI mask
#' @export

ndvi.mask <- function(hy.file, metadata.path, reflectance.path, wavelength.path, red.nm, nir.nm, ndvi.threshold){
 
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
  
  
  # lets find the index for the closest band to NIR
  red.wl <- red.nm
  red.index <- which(abs(wavelengths - red.wl) == min(abs(wavelengths - red.wl)))
  
  # lets find the index for the closest band to Red
  nir.wl <- nir.nm
  nir.index <- which(abs(wavelengths - nir.wl) == min(abs(wavelengths - nir.wl)))
  
  # lets read in the NIR band and clean it up like we need before
  nir.array <- h5read(file = hy.file,
                      name = reflectance.path,
                      index = list(nir.index, 1:n.cols, 1:n.rows))
  
  nir.matrix <- nir.array[1,,]
  nir.matrix[nir.matrix == data.ignore.val] <- NA
  nir.matrix <- nir.matrix / scale.fact.val
  
  # lets read in the RED band and clean it up
  red.array <- h5read(file = hy.file,
                      name = reflectance.path,
                      index = list(red.index, 1:n.cols, 1:n.rows))
  
  red.matrix <- red.array[1,,]
  red.matrix[red.matrix == data.ignore.val] <- NA
  red.matrix <- red.matrix / scale.fact.val
  
  # lets make a NDVI array
  ndvi.array <- (nir.matrix - red.matrix) / (nir.matrix + red.matrix)
  
  # lets make a mask where NDVI < the given threshold = 0
  ndvi.array[ndvi.array < ndvi.threshold] <- 0
  ndvi.array[ndvi.array >= ndvi.threshold] <- 1
  
  # return the array
  return(ndvi.array)
}