## example marker df
marker <- data.frame(
    markerID = c("marker1", "marker2"),
    rtime = c(20, 80),
    mobility = c(0, 2000)
)

test_that("Transformation of numeric vectors works", {
    rtime <- c(10, 20, 30, 40, 50, 60, 70, 80, 90, 100)
    transf <- mobilityTransform(x = rtime, marker = marker)

    expect_true(is.vector(transf))
    expect_true(is.numeric(transf))
    expect_equal(length(transf), length(rtime))
    expect_equal(transf[1], -2666.667, tolerance = 1e-06)
    expect_equal(transf[2], 0)
    expect_equal(transf[3], 888.8889, tolerance = 1e-06)
    expect_equal(MobilityTransformR:::.transformNumeric(rtime, marker,
        tR = 0,
        U = 30, L = 90
    ), transf)
    expect_error(
        mobilityTransform(rtime, marker[1, ]),
        "'U' and 'L' are expected to be of length 1."
    )
    expect_error(
        mobilityTransform(c(1, 2, 3, "a"), marker),
        "'x' needs to be of class 'numeric', 'Spectra' or 'OnDistMSnExp' but not class 'character'"
    )
    expect_error(
        mobilityTransform(rtime),
        "Missing data.frame 'marker' with marker information"
    )
    expect_error(
        mobilityTransform(rtime, "A"),
        "'marker' missing column 'rtime', 'mobility' or both"
    )
})

test_that("Transformation of Spectra works", {
    spectra_data <- Spectra::Spectra(system.file("extdata/CEMS_10ppm.mzML",
        package = "MobilityTransformR"
    ))

    transf <- mobilityTransform(x = spectra_data, marker = marker)

    expect_true(is(transf, "Spectra"))
    expect_true(is.numeric(transf$rtime))
    expect_equal(length(transf$rtime), length(spectra_data$rtime))
    expect_equal(transf$rtime[1], 2533.447, tolerance = 1e-06)
    expect_equal(transf$rtime[10], 2534.953, tolerance = 1e-06)
    expect_equal(transf$rtime[100], 2548.345, tolerance = 1e-06)
    expect_equal(MobilityTransformR:::.transformSpectra(spectra_data, marker,
        tR = 0,
        U = 30, L = 90
    ), transf)
    expect_error(
        mobilityTransform(spectra_data, marker[1, ]),
        "'U' and 'L' are expected to be of length 1."
    )
    expect_error(
        mobilityTransform("A", marker),
        "'x' needs to be of class 'numeric', 'Spectra' or 'OnDistMSnExp' but not class 'character'"
    )
    expect_error(
        mobilityTransform(spectra_data),
        "Missing data.frame 'marker' with marker information"
    )
    expect_error(
        mobilityTransform(spectra_data, "A"),
        "'marker' missing column 'rtime', 'mobility' or both"
    )

    ## test if MSnbase/ xcms functions work
    fl <- tempfile()
    Spectra::export(transf, Spectra::MsBackendMzR(), file = fl)
    ## load again as OnDiskMSnExp
    transf_load <- readMSData(files = fl, mode = "onDisk")
    mz_test <- c(152.071154 - 0.001, 152.071154 + 0.001)
    EIE_test <- transf_load |> MSnbase::filterMz(mz = mz_test)
    xdata <- xcms::findChromPeaks(EIE_test,
        param = xcms::MatchedFilterParam(binSize = 1)
    )

    c_peaks <- as.data.frame(chromPeaks(xdata))
    expect_true(is(xdata, "XCMSnExp"))
    expect_true(is.data.frame(c_peaks))
    expect_equal(nrow(c_peaks), 1)
    expect_equal(c_peaks$rt, 2603.273, tolerance = 1e-06)
})


test_that("Transformation of OnDiskExp works", {
    fl <- system.file("extdata/CEMS_10ppm.mzML",
        package = "MobilityTransformR"
    )
    raw_data <- MSnbase::readMSData(
        files = fl,
        mode = "onDisk"
    )

    marker <- data.frame(
        markerID = c("marker1", "marker2"),
        rtime = c(20, 80),
        mobility = c(0, 2000),
        fileIdx = c(1, 1)
    )

    transf <- mobilityTransform(x = raw_data, marker = marker)

    expect_true(is(transf, "OnDiskMSnExp"))
    expect_true(is.numeric(rtime(transf)))
    expect_equal(length(rtime(transf)), length(rtime(raw_data)))
    expect_equal(rtime(transf)[[1]], 2533.447, tolerance = 1e-06)
    expect_equal(rtime(transf)[[10]], 2534.953, tolerance = 1e-06)
    expect_equal(rtime(transf)[[100]], 2548.345, tolerance = 1e-06)
    expect_equal(MobilityTransformR:::.transformOnDiskMSnExp(raw_data, marker,
        tR = 0, U = 30, L = 90
    ), transf)
    expect_error(
        mobilityTransform(raw_data, marker[1, ]),
        "'U' and 'L' are expected to be of length 1."
    )
    expect_error(
        mobilityTransform("A", marker),
        "'x' needs to be of class 'numeric', 'Spectra' or 'OnDistMSnExp' but not class 'character'"
    )
    expect_error(
        mobilityTransform(raw_data),
        "Missing data.frame 'marker' with marker information"
    )
    expect_error(
        mobilityTransform(raw_data, "A"),
        "'marker' missing column 'rtime', 'mobility' or both"
    )
    expect_error(
        mobilityTransform(raw_data, marker[, -4]),
        "'marker' missing column 'fileIdx'"
    )

    ## test if MSnbase/ xcms functions work
    ### test export
    fl <- tempfile()
    MSnbase::writeMSData(transf, file = fl, copy = FALSE)

    ## load again as OnDiskMSnExp
    transf_load <- readMSData(files = fl, mode = "onDisk")
    mz_test <- c(152.071154 - 0.001, 152.071154 + 0.001)
    EIE_test <- transf_load |> MSnbase::filterMz(mz = mz_test)
    xdata <- xcms::findChromPeaks(EIE_test,
        param = xcms::MatchedFilterParam(binSize = 1)
    )

    c_peaks <- as.data.frame(chromPeaks(xdata))
    expect_true(is(xdata, "XCMSnExp"))
    expect_true(is.data.frame(c_peaks))
    expect_equal(nrow(c_peaks), 1)
    expect_equal(c_peaks$rt, 2603.273, tolerance = 1e-06)
})
