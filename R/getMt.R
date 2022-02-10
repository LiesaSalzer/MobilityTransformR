#' @name getMtime
#'
#' @aliases getMtime
#'
#' @title Get the migration time of a known compound 
#'
#' @description
#' The function `getMtime` searches the migration time (mt) of a known compound 
#' within a specified mz and mt range (`mz` and `mt`). 
#'
#' @param x
#' `OnDiskMSnExp` object, containing CE-MS(/MS) data from a single file (e.g.
#' filtered by `filterFile()`). 
#' The migration time (mt) is provided in sec. The CE-MS data contain known 
#' compounds `mz`, of which the exact migration time is going to be determined. 
#'
#' @param mz
#' `numeric`, with upper and lower limit, (length(`mz`) = 2), representing the 
#' mz range of the compound of which the migration time is going to be 
#' determined. The range should be as narrow as possible but depends on the mass 
#' accuracy of the mass spectrometer that has been used to acquire the data. 
#'
#' @param mt
#' `numeric`, with upper and lower limit, (length(`mt`) = 2), limiting the 
#' migration tame range of the compound.
#' Use a narrow `mt` range, to to avoid that other components with the same mz 
#' and different mt are being picked. 
#' 
#' @param param
#' `method`, from xcms that defined how peaks will be picked. The default is 
#' `MatchedFilterParam(binSize = 1, snthresh = 100)`
#'
#' @details
#' `getMtime` uses CE-MS data stored in `OnDiskMSnExp` objects to search for the 
#' migration time of selected compounds as for example EOF markers in order to 
#' perform effective mobility scale transformation. 
#' The `OnDiskMSnExp` object is filtered using the defined mz-range `mz`, and 
#' the mt-range `mt`, where the compound is expected to migrate. 
#' The migration time is determined by applying the peak picking algorithm from
#' xcms. 
#'
#' @return
#' `data.frame` with two columns "rtime" storing the migration time in sec, 
#' and "fileIdx" storing the file index and the number of rows corresponding to 
#' the number of input files. 
#'
#' @author Liesa Salzer 
#'
#'  
#'@importFrom  xcms findChromPeaks 
#'@importFrom  xcms chromPeaks
#'
#' @examples 
#' fl <- system.file("extdata/CEMS_10ppm.mzML", 
#' package = "MobilityTransformR")
#' raw_data <- readMSData(files = fl,
#'                        mode = "onDisk")
#'# [M+H]+ of paracetamol: mz = 152.071154
#' mz_paracetamol <- c(152.071154 - 0.005, 152.071154 + 0.005)
#' mt_paracetamol <- c(600, 900)                      
#' getMtime(raw_data, mz = mz_paracetamol, mt = mt_paracetamol)                       
#'                        
#'@export
getMtime <- function(x, mz = numeric(), mt = numeric(), 
                     param = MatchedFilterParam(binSize = 1, snthresh = 50),
                     ...) {
## sanity checks
if (!is(x, "OnDiskMSnExp")) 
  stop("'x' is not of class 'OnDiskMSnExp'!")

if (missing(mz) | missing(mt)) 
  stop("Arguments 'mz' and 'mt' are required!")


  rt_df <- data.frame(rtime = integer(),
                      fileIdx = integer())
  
  for (i in unique(fromFile(x))) {
    x_i <- filterFile(x, i)
    
    df <- findChromPeaks(filterMz(x_i, mz = mz) |> filterRt(rt = mt), 
      param = param
    ) |> chromPeaks() |> 
      as.data.frame()
    
    
    if (length(df$rt) != 1) 
      stop(length(df$rt), " peaks have been found in file ", i,
           ", align input parameters")
    
    else 
      rt_df <- rbind(rt_df, data.frame(rtime = df$rt,
                                      fileIdx = i))

    
  }
  
  return(rt_df)


}
