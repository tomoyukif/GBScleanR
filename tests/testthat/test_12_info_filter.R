# sample.gds ships with MQ and QD INFO fields used by setInfoFilter().
gds_fn <- system.file("extdata", "sample.gds", package = "GBScleanR")

test_that("getInfo returns INFO annotations", {
    gds <- loadGDS(gds_fn, verbose = FALSE)
    mq <- getInfo(gds, "MQ")
    qd <- getInfo(gds, "QD")
    expect_equal(length(mq), nmar(gds))
    expect_equal(length(qd), nmar(gds))
    expect_true(is.numeric(mq))
    expect_true(is.numeric(qd))

    valid_snp <- sample(c(TRUE, FALSE), nmar(gds), TRUE)
    validMar(gds) <- valid_snp
    expect_equal(getInfo(gds, "MQ", valid = TRUE), mq[valid_snp],
                 ignore_attr = TRUE)
    closeGDS(gds, verbose = FALSE)
})

test_that("setInfoFilter validates arguments", {
    gds <- loadGDS(gds_fn, verbose = FALSE)
    expect_error(setInfoFilter(gds, mq = -1), "mq should be")
    expect_error(setInfoFilter(gds, fs = -1), "fs should be")
    expect_error(setInfoFilter(gds, qd = -1), "qd should be")
    expect_error(setInfoFilter(gds, mqranksum = 1), "mqranksum should be")
    expect_error(setInfoFilter(gds, mqranksum = c(1, 0)),
                 "mqranksum\\[1\\] should be smaller")
    closeGDS(gds, verbose = FALSE)
})

test_that("setInfoFilter filters markers by MQ and QD", {
    gds <- loadGDS(gds_fn, verbose = FALSE)
    mq <- getInfo(gds, "MQ")
    qd <- getInfo(gds, "QD")

    gds <- setInfoFilter(gds, mq = 30)
    expect_equal(validMar(gds), mq >= 30, ignore_attr = TRUE)

    gds <- resetMarFilter(gds)
    gds <- setInfoFilter(gds, qd = 30)
    expect_equal(validMar(gds), qd >= 30, ignore_attr = TRUE)

    gds <- resetMarFilter(gds)
    gds <- setInfoFilter(gds, mq = 25, qd = 25)
    expect_equal(validMar(gds), mq >= 25 & qd >= 25, ignore_attr = TRUE)
    closeGDS(gds, verbose = FALSE)
})
