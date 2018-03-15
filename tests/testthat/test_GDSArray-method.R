test_that("gdsfile getter and setter works", {
    file <- system.file(package="SNPRelate", "extdata", "hapmap_geno.gds")
    gds <- GDSArray(file)
    expect_equal(gdsfile(gds), file)
    file1 <- SeqArray::seqExampleFileName("gds")
    gdsfile(gds) <- file1
    expect_true(validObject(gds))
})

