test_that("get GDS file format works", {
    .get_gds_fileFormat <- GDSArray:::.get_gds_fileFormat
    file <- system.file(package="SNPRelate", "extdata", "hapmap_geno.gds")
    ff <- .get_gds_fileFormat(file)
    expect_true(is.character(ff))
    expect_true(ff %in% c("SNP_ARRAY", "SEQ_ARRAY"))

    file1 <- system.file(package="SeqArray", "extdata", "CEU_Exon.gds")
    ff1 <- .get_gds_fileFormat(file1)
    expect_true(is.character(ff1))
    expect_true(ff1 %in% c("SNP_ARRAY", "SEQ_ARRAY"))
})

test_that("get gds all nodes works", {
    file <- system.file(package="SeqArray", "extdata", "CEU_Exon.gds")
    ## allnodes
    allnodes <- gdsnodes(file)
    expect_true(is(allnodes, "character"))
})

test_that("extracting gds node dimension orientation works", {
    .read_gdsnode_sampleInCol <- GDSArray:::.read_gdsnode_sampleInCol

    file <- system.file(package="SNPRelate", "extdata", "hapmap_geno.gds")
    val <- .read_gdsnode_sampleInCol(file, "genotype")
    expect_false(val)
    
    file1 <- system.file(package="SeqArray", "extdata", "CEU_Exon.gds")
    expect_false(.read_gdsnode_sampleInCol(file1, "genotype"))
    expect_false(.read_gdsnode_sampleInCol(file1, "annotation/format/DP"))
})


test_that("extracting gds node dimnames works", {
    .get_gdsnode_dim <- GDSArray:::.get_gdsnode_dim
    .get_gdsnode_dimnames <- GDSArray:::.get_gdsnode_dimnames

    file <- system.file(package="SeqArray", "extdata", "CEU_Exon.gds")
    dimnames <- .get_gdsnode_dimnames(file, "genotype")
    expect_true(is(dimnames, "list"))
    expect_true(all( vapply(dimnames, is, logical(1), "character") ))

    dims <- .get_gdsnode_dim(file, "genotype")
    expect_true(identical(unname(lengths(dimnames)), dims))
})
