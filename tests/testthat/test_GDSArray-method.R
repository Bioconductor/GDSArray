test_that("gdsfile getter and setter works", {
    file <- system.file(package="SNPRelate", "extdata", "hapmap_geno.gds")
    ga <- GDSArray(file)
    expect_equal(gdsfile(ga), file)
    file1 <- SeqArray::seqExampleFileName("gds")
    gdsfile(ga) <- file1
    expect_true(validObject(ga))
})

