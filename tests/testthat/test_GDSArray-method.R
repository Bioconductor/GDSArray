test_that("gdsfile getter and setter works", {
    file <- gdsExampleFileName("snpgds")
    ga <- GDSArray(file, "genotype")
    expect_equal(gdsfile(ga), file)
    file1 <- gdsExampleFileName("seqgds")
    gdsfile(ga) <- file1
    expect_true(validObject(ga))
})

