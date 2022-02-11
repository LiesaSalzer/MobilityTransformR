## load toy example data set
fl <- system.file("extdata/CEMS_10ppm.mzML",
                  package = "MobilityTransformR")
raw_data <- MSnbase::readMSData(files = fl,
                       mode = "onDisk")
# [M+H]+ of paracetamol: mz = 152.071154
mz_paracetamol <- c(152.065154, 152.076154)
mt_paracetamol <- c(600, 900)


test_that("Getting migration time works", {
  mt <- getMtime(raw_data, mz = mz_paracetamol, mt = mt_paracetamol)
  expect_true(is.data.frame(mt))
  expect_equal(colnames(mt), c("rtime", "fileIdx"))
  expect_equal(dim(mt), c(1, 2))
  expect_true(is.numeric(mt$rtime))
  expect_true(is.numeric(mt$fileIdx))
  expect_equal(mt[1,1], 840.796, tolerance = 1e-06)
  expect_equal(mt[1,2], 1)
  
  expect_error(getMtime(NA, mz = mz_paracetamol, mt = mt_paracetamol),
               "'x' is not of class 'OnDiskMSnExp'!")
  expect_error(getMtime(raw_data, mz = mz_paracetamol),
               "Arguments 'mz' and 'mt' are required!")
  expect_error(getMtime(raw_data, mt = mt_paracetamol),
               "Arguments 'mz' and 'mt' are required!")
  expect_error(getMtime(raw_data, mz = c(150,155), mt = mt_paracetamol),
               "3 peaks have been found in file 1, align input parameters")
  expect_warning(expect_error(getMtime(raw_data, mz = c(155.0001,155.0002), mt = mt_paracetamol),
               "0 peaks have been found in file 1, align input parameters"))
  
})

test_that("Getting migration time with CentWaveParam works", {
  mt_cwp <- getMtime(raw_data, mz = c(152,152.2), mt = mt_paracetamol, 
                     param = xcms::CentWaveParam())
  
  expect_true(is.data.frame(mt_cwp))
  expect_equal(colnames(mt_cwp), c("rtime", "fileIdx"))
  expect_equal(dim(mt_cwp), c(1, 2))
  expect_true(is.numeric(mt_cwp$rtime))
  expect_true(is.numeric(mt_cwp$fileIdx))
  expect_equal(mt_cwp[1,1], 841.305, tolerance = 1e-06)
  expect_equal(mt_cwp[1,2], 1)
  
  expect_error(getMtime(NA, mz = mz_paracetamol, mt = mt_paracetamol),
               "'x' is not of class 'OnDiskMSnExp'!")
  expect_error(getMtime(raw_data, mz = mz_paracetamol),
               "Arguments 'mz' and 'mt' are required!")
  expect_error(getMtime(raw_data, mt = mt_paracetamol),
               "Arguments 'mz' and 'mt' are required!")
  expect_error(getMtime(raw_data, mz = c(150,155), mt = mt_paracetamol),
               "3 peaks have been found in file 1, align input parameters")
  expect_warning(expect_error(getMtime(raw_data, mz = c(155.0001,155.0002), 
                                       mt = mt_paracetamol),
               "0 peaks have been found in file 1, align input parameters"))
})
