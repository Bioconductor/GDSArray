test_that("get GDS file format works", {
    .get_gdsdata_fileFormat <- GDSArray:::.get_gdsdata_fileFormat
    file <- system.file(package="SNPRelate", "extdata", "hapmap_geno.gds")
    ff <- .get_gdsdata_fileFormat(file)
    expect_true(is.character(ff))
    expect_true(ff %in% c("SNP_ARRAY", "SEQ_ARRAY"))

    file1 <- system.file(package="SeqArray", "extdata", "CEU_Exon.gds")
    ff1 <- .get_gdsdata_fileFormat(file1)
    expect_true(is.character(ff1))
    expect_true(ff1 %in% c("SNP_ARRAY", "SEQ_ARRAY"))
})

test_that("get gds all nodes, array nodes, and node dimension works", {
    .get_gdsdata_arrayNodes <- GDSArray:::.get_gdsdata_arrayNodes
    .get_gdsdata_dim <- GDSArray:::.get_gdsdata_dim
    
    file <- system.file(package="SeqArray", "extdata", "CEU_Exon.gds")
    ## allnodes
    allnodes <- gdsNodes(file)
    expect_true(is(allnodes, "character"))

    ## dims
    dims <- lapply(arraynodes, function(x) .get_gdsdata_dim(file, x))
    expect_true(validObject(dims))

    ## array nodes
    ndims <- lengths(dims)
    expect_true(all(ndims > 1))

    
})

test_that("extracting gds node dimension orientation works", {
    .read_gdsdata_sampleInCol <- GDSArray:::.read_gdsdata_sampleInCol

    file <- system.file(package="SNPRelate", "extdata", "hapmap_geno.gds")
    val <- .read_gdsdata_sampleInCol(file, "genotype")
    expect_false(val)
    
    file1 <- system.file(package="SeqArray", "extdata", "CEU_Exon.gds")
    expect_false(.read_gdsdata_sampleInCol(file1, "genotype/data"))
    expect_false(.read_gdsdata_sampleInCol(file1, "annotation/format/DP/data"))
    expect_false(.read <- gdsdata <- sampleInCol(file1, "phase/data"))
})


test_that("extracting gds node dimnames works", {
    .get_gdsdata_dim <- GDSArray:::.get_gdsdata_dim
    .get_gdsdata_dimnames <- GDSArray:::.get_gdsdata_dimnames

    file <- system.file(package="SeqArray", "extdata", "CEU_Exon.gds")
    dimnames <- .get_gdsdata_dimnames(file, "genotype/data")
    expect_true(is(dimnames, "list"))

    dims <- .get_gdsdata_dim(file, "genotype/data")
    expect_true(identical(unname(lengths(dimnames)), dims))
})
