.get_gds_fileFormat <- function(file)
{
    f <- openfn.gds(file)
    on.exit(closefn.gds(f))
    ff <- get.attr.gdsn(f$root)$FileFormat
    ff
}

## .get_gdsnode_isarray <- function(file, node)
## {
##     f <- openfn.gds(file)
##     on.exit(closefn.gds(f))
    
##     isarray <- objdesp.gdsn(index.gdsn(f, node))$is.array
##     return(isarray)
## }

## ## array data with >1 dimensions to pass into assays(se)
## .get_gdsnode_non1D_array <- function(file)
## {
##     names.gdsn <- gdsnodes(file)
##     isarray <- vapply(names.gdsn,
##                       function(x) .get_gdsnode_isarray(file, x),
##                       logical(1))
##     dims <- lapply(names.gdsn, function(x) .get_gdsnode_dim(file, x))
##     names.gdsn[
##         isarray & lengths(dims) > 1 & 
##         ! vapply(dims, function(x) any(x == 0L), logical(1)) &
##         !grepl("~", names.gdsn)
##         ## what's the pattern with genotype/~data, phase/~data,
##         ## annotation/format/DP/~data
##     ]
## }

.read_gdsnode_sampleInCol <- function(file, node)
{
    ff <- .get_gds_fileFormat(file)
    dims <- .get_gdsnode_dim(file, node)
    stopifnot(length(dims) > 1)
    
    if (ff == "SNP_ARRAY") {
        f <- openfn.gds(file)
        on.exit(closefn.gds(f))
        
        rd <- names(get.attr.gdsn(index.gdsn(f, node)))
        if ("snp.order" %in% rd) sampleInCol <- TRUE   ## snpfirstdim (in row)
        if ("sample.order" %in% rd) sampleInCol <- FALSE
    } else if (ff == "SEQ_ARRAY") {
        f <- seqOpen(file)
        on.exit(seqClose(f))
        seqSumm <- seqSummary(f, verbose=FALSE)
        dimSumm <- c(
            ploidy = seqSumm$ploidy,
            sample = seqSumm$num.sample,
            variant = seqSumm$num.variant)
        ind <- match(dimSumm[c("variant", "sample")], dims)
        if (ind[1] < ind[2]) {
            sampleInCol <- TRUE   ## rewrite for general cases: format/DP.
        } else {
            sampleInCol <- FALSE
        }
    }
    sampleInCol
}

.get_gdsnode_type <- function(file, node)
{
    f <- openfn.gds(file)
    on.exit(closefn.gds(f))
    
    type <- objdesp.gdsn(index.gdsn(f, node))$type
    type_levels <- c(
        "Label", "Folder", "VFolder", "Raw", "Integer",
        "Factor", "Logical", "Real", "String", "Unknown")
    if (!type %in% type_levels)
        stop(
            wmsg("The type of GDS nodes should be one of:\n",
                 paste0(type_levels, collapse=", ")))
    type <- as.character(type)
    type
}

.get_gdsnode_dim <- function(file, node)
{
    ff <- .get_gds_fileFormat(file)
    dimSumm <- lengths(.get_gds_dimnames(file))
    f <- openfn.gds(file)
    on.exit(closefn.gds(f))
    if (ff == "SEQ_ARRAY" && grepl("genotype$|DP$", node)) {
        node <- paste0(node, "/data")
    }
    dim <- objdesp.gdsn(index.gdsn(f, node))$dim
    ## }
    if (is.null(dim))
        return(NULL)
    ind <- match(dim, dimSumm)  ## make sure the dim aligns with
                                ## dimnames lengths
    if (length(ind) ==1 && is.na(ind)) {  ## in cases like "annotation/info/AA" length
                       ## 1328 not 1348.
        if (grep("^annotation", node))
            dim <- unname(dimSumm["variant.id"])
    }
    if (!is.integer(dim)) {
        if (any(dim > .Machine$integer.max)) {
            dim_in1string <- paste0(dim, collapse=" x ")
            stop(wmsg(
                "The dimensions of GDS dataset '", file, "' are: ",
                dim_in1string, "\n\nThe GDSArray package only ",
                "supports datasets with all dimensions <= 2^31-1",
                " (this is ", .Machine$integer.max, ") at the moment."))
        }
    }
    dim <- as.integer(dim)
    dim
}

.get_gdsnode_dimnames <- function(file, node)
{
    dims <- .get_gdsnode_dim(file, node)  ## for "SEQ_ARRAY", do not use objdesp for dim.
    dimnames <- .get_gds_dimnames(file)
    dimSumm <- lengths(dimnames)
    if(is.null(dims) | any(dims == 0))
        return(NULL)
    ind <- match(dims, dimSumm)
    dimnames <- dimnames[ind]
    if (!is.list(dimnames)) {
        stop(wmsg(
            "The dimnames of GDS dataset '", file,
            "' should be a list!"))
    }
    dimnames
}

## dimnames for GDSArraySeed must be character!!!
## http://bioconductor.org/packages/release/bioc/vignettes/DelayedArray/inst/doc/02-Implementing_a_backend.html#implementing-optimized-backend-specific-methods
.get_gds_dimnames <- function(file)  
{
    ff <- .get_gds_fileFormat(file)
    
    if (ff == "SNP_ARRAY") {
        file.summ <- snpgdsSummary(file, show=FALSE)
        dimnames <- lapply(file.summ[c("sample.id", "snp.id")], as.character)
        ## f <- openfn.gds(file)
        ## on.exit(closefn.gds(f))
        
        ## sample.id <- read.gdsn(index.gdsn(f, "sample.id"))
        ## snp.id <- read.gdsn(index.gdsn(f, "snp.id"))
        ## rd <- names(get.attr.gdsn(index.gdsn(f, node))) ## ?
        ## if ("snp.order" %in% rd) {
        ##     dimnames <-
        ##         list(snp.id = as.character(snp.id), sample.id = sample.id)
        ## } else {
        ##     dimnames <-
        ##         list(sample.id = sample.id, snp.id = as.character(snp.id))
        ## }
        ## dimSumm <- lengths(dimnames)
    } else if (ff == "SEQ_ARRAY") {
        f <- seqOpen(file)
        on.exit(seqClose(f))
        sample.id <- read.gdsn(index.gdsn(f, "sample.id"))
        variant.id <- read.gdsn(index.gdsn(f, "variant.id"))
        file.summ <- seqSummary(f, verbose=FALSE)
        ## dimSumm <- c(
        ##     ploidy = file.summ$ploidy,
        ##     sample = file.summ$num.sample,
        ##     variant = file.summ$num.variant)
        ## stopifnot(length(variant.id) == dimSumm["variant"])
        ## stopifnot(length(sample.id) == dimSumm["sample"])
        dimnames <- list(
            ploidy.id = as.character(seq_len(file.summ$ploidy)),
            sample.id = as.character(sample.id),
            variant.id = as.character(variant.id))
    }
    return(dimnames)
}
## FIXME: if dimnames is shorter than the corresponding dimensions,
## use NULL to extend.??

.get_gdsnode_first_val <- function(file, node)
{
    ff <- .get_gds_fileFormat(file)
    dims <- .get_gdsnode_dim(file, node)
    if (is.null(dims) | any(dims == 0L))
        return(NULL)
    if (ff == "SEQ_ARRAY") {
        f <- seqOpen(file)
        on.exit(seqClose(f))
        seqSetFilter(f, variant.sel = 1, sample.sel = 1, verbose = FALSE)
        first_val <- seqGetData(f, node)
        if(is.list(first_val)) {  ## for "annotation/format/DP" only
            first_val <- first_val$data[1,1]
        } else if (length(dim(first_val)) == 3) { ## for "genotype" only
            first_val <- first_val[1,,]
        }
        seqResetFilter(f, verbose = FALSE)
    } else {
        f <- openfn.gds(file)
        on.exit(closefn.gds(f))
        first_val <- readex.gdsn(
            index.gdsn(f, node), sel=as.list(rep(1, length(dims))))
    }
    first_val
}


example <- function(name="GDSArray")
{
    file <- SeqArray::seqExampleFileName("gds")
    gds <- GDSArray(file, "annotation/format/DP/data")
    gds
}
