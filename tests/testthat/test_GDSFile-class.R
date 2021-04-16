file <- gdsExampleFileName("seqgds")
file1 <- gdsExampleFileName("snpgds")

test_that("GDSFile constructor works", {    
    gf <- GDSFile(file)
    expect_s4_class(gf, "GDSFile")
    expect_true(validObject(gf))
    
})

test_that("gdsfile getter and setter works", {
    gf <- GDSFile(file)
    expect_equal(gdsfile(gf), file)
    
    gdsfile(gf) <- file1
    expect_true(validObject(gf))
})

test_that("$ completion works", {
    gf <- GDSFile(file)
    .DollarNames.GDSFile <- GDSArray:::.DollarNames.GDSFile
    expect_true("annotation" %in% .DollarNames.GDSFile(gf, "anno"))
    expect_true("info" %in% .DollarNames.GDSFile(gf$annotation, "in"))
})

test_that("$ works", {
    gf <- GDSFile(file)
    expect_true(validObject(gf$sample.id))
    expect_s4_class(gf$annotation, "GDSFile")
    expect_s4_class(gf$annotation$info$AC, "GDSArray")
})

test_that("gdsnodes works", {
    gf <- GDSFile(file)
    expect_true(all(gdsnodes(gf) %in% gdsnodes(file)))
})
