#' Convert a hdf5 file containing hyperspectral imagery into a brightness mask for use in later corrections.
#'
#' This function reads in a hdf5 file and create a brightness matrix to remove shaded pixels for use in the 
#' BRDF correction that is applied in later steps.
#' 
#' For theory behind this function please see:
#' 
#' Clark, M.L., Roberts, D.A., and Clark, D.B., 2005. Hyperspectral discrimination of tropical rain forest
#' tree species at crown scales. Remote Sensing of Environment, 96: 375-398.
#' 
#' Gougeon, F.A., 1995. Comparison of possible multispectral classification schemes for tree crowns individually
#' delineated on high spatial resolution MEIS images. Canadian Journal of Remote Sensing, 21(1): 1-9.
#'
#' @param hy.file hdf5 file containing hyperspectral imagery and associated metadata
#' @param metadata.path hdf5 path to reflectance metadata
#' @param reflectance.path hdf5 path to reflectance data
#' @param wavelength.path hdf5 path to wavelength metadata
#' @return A matrix containing the brightness mask
#' @export

brightness.mask <- function(hy.file, metadata.path, reflectance.path, wavelength.path){
  
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
  
  # lets find the index for the closest band to 800nm for NIR
  wl <- 800
  nir.800.index <- which(abs(wavelengths - wl) == min(abs(wavelengths - wl)))
  
  # lets read in a single band from one of the hdf5 files
  nir.800.array <- h5read(file = hy.file,
                          name = reflectance.path,
                          index = list(nir.800.index, 1:n.cols, 1:n.rows))
  
  # lets convert to a raster
  nir.800.matrix <- nir.800.array[1,,]
  nir.800.matrix[nir.800.matrix == data.ignore.val] <- NA
  nir.800.mat <- nir.800.matrix / scale.fact.val
  
  # lets calculate the mean reflectance value for this matrix
  nir.800.mean <- mean(nir.800.mat, na.rm = TRUE)
  
  # lets make a mask where reflectance < the mean = 0
  nir.800.mat[nir.800.mat < nir.800.mean] <- 0
  nir.800.mat[nir.800.mat >= nir.800.mean] <- 1
  
  # return the array
  return(nir.800.mat)
}