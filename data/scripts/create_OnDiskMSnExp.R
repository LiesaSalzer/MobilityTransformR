#' Creating the OnDiskMSnExp from .mzML raw data.
library(MSnbase)

fl <- list.files(system.file("extdata/", 
                             package = "MobilityTransformR"), 
                 pattern = ".mzML", full.names = T)

# Load mzXML data with MSnBase
raw_data_all <- readMSData(files = fl, mode = "onDisk")
raw_data <- filterFile(raw_data_all, 1)

save(raw_data, file = "data/CEMS_OnDisk.RData")
