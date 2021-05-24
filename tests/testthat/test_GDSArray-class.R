file <- gdsExampleFileName("snpgds")

test_that("GDSArraySeed constructor works", {
    f <- openfn.gds(file)
    seed <- GDSArray:::GDSArraySeed(f, "genotype")
    expect_s4_class(seed, "GDSArraySeed")
    expect_true(validObject(seed))
    ## expect_equal(dim(seed), c(9088L, 279L))
    expect_equal(dim(seed), c(279L, 9088L))
    expect_equal(class(dimnames(seed)), "list")
    ## expect_equal(lengths(unname(dimnames(seed))), dim(seed))
    expect_equal(dimnames(seed), list(NULL, NULL))
    closefn.gds(f)
})


test_that("GDSArray constructor works", {
    ga <- GDSArray(file, "genotype")
    expect_s4_class(ga, "GDSArray")
    expect_true(validObject(ga))
    ## expect_equal(dim(ga), c(9088L, 279L))
    expect_equal(dim(ga), c(279L, 9088L))
    seed <- seed(ga)
    expect_s4_class(seed, "GDSArraySeed")
    expect_true(validObject(seed))
})

