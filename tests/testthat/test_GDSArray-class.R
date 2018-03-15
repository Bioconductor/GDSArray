test_that("GDSArraySeed constructor works", {
    file <- system.file(package="SNPRelate", "extdata", "hapmap_geno.gds")
    seed <- GDSArraySeed(file, "genotype")
    expect_s4_class(seed, "GDSArraySeed")
    expect_true(validObject(seed))
    expect_equal(dim(seed), c(9088L, 279L))
    expect_equal(class(dimnames(seed)), "list")
    expect_equal(lengths(unname(dimnames(seed))), dim(seed))
})


test_that("GDSArray constructor works", {
    file <- system.file(package="SNPRelate", "extdata", "hapmap_geno.gds")
    gds <- GDSArray(file)
    expect_s4_class(gds, "GDSArray")
    expect_true(validObject(gds))
    expect_equal(dim(gds), c(9088L, 279L))
    seed <- seed(gds)
    expect_s4_class(seed, "GDSArraySeed")
    expect_true(validObject(seed))
})

