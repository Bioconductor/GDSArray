test_that("GDSFile constructor works", {
    file <- SeqArray::seqExampleFileName("gds")
    gf <- GDSFile(file)
    expect_s4_class(gf, "GDSFile")
    expect_true(validObject(gf))
    
})

test_that("gdsfile getter and setter works", {
    file <- SeqArray::seqExampleFileName("gds")
    gf <- GDSFile(file)
    expect_equal(gdsfile(gf), file)
    
    file1 <- SNPRelate::snpgdsExampleFileName()
    gdsfile(gf) <- file1
    expect_true(validObject(gf))
})

test_that("$ completion works", {
    file <- SeqArray::seqExampleFileName("gds")
    gf <- GDSFile(file)
    .DollarNames.GDSFile <- GDSArray:::.DollarNames.GDSFile
    expect_true("annotation" %in% .DollarNames.GDSFile(gf, "anno"))
    expect_true("info" %in% .DollarNames.GDSFile(gf$annotation, "in"))
})

test_that("$ works", {
    file <- SeqArray::seqExampleFileName("gds")
    gf <- GDSFile(file)
    expect_true(validObject(gf$sample.id))
    expect_s4_class(gf$annotation, "GDSFile")
    expect_s4_class(gf$annotation$info$AC, "GDSArray")
})

test_that("gdsnodes works", {
    file <- SeqArray::seqExampleFileName("gds")
    gf <- GDSFile(file)
    expect_identical(gdsnodes(file), gdsnodes(gf))
})
