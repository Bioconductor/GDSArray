test_that("GDSArraySeed constructor works", {
    file <- system.file(package="SNPRelate", "extdata", "hapmap_geno.gds")
    seed <- GDSArray:::GDSArraySeed(file, "genotype")
    expect_s4_class(seed, "GDSArraySeed")
    expect_true(validObject(seed))
    expect_equal(dim(seed), c(9088L, 279L))
    expect_equal(class(dimnames(seed)), "list")
    expect_equal(lengths(unname(dimnames(seed))), dim(seed))
})


test_that("GDSArray constructor works", {
    file <- system.file(package="SNPRelate", "extdata", "hapmap_geno.gds")
    ga <- GDSArray(file)
    expect_s4_class(ga, "GDSArray")
    expect_true(validObject(ga))
    expect_equal(dim(ga), c(9088L, 279L))
    seed <- seed(ga)
    expect_s4_class(seed, "GDSArraySeed")
    expect_true(validObject(seed))
})

