.get_gdsdata_fileFormat <- function(file)
{
    f <- openfn.gds(file)
    on.exit(closefn.gds(f))
    ff <- get.attr.gdsn(f$root)$FileFormat
    ff
}

## array data with >1 dimensions to pass into assays(se)
.get_gdsdata_arrayNodes <- function(file)
{
    names.gdsn <- gdsNodes(file)
    f <- openfn.gds(file)
    on.exit(closefn.gds(f))
    
    isarray <- vapply(names.gdsn, function(x) {
        objdesp.gdsn(index.gdsn(f, x))$is.array
    }, logical(1))
    dims <- lapply(names.gdsn, function(x)
        objdesp.gdsn(index.gdsn(f, x))$dim)
    ## names(dims) <- all.gdsn
    names.gdsn[
        isarray & lengths(dims) > 1 & 
        ! unlist(lapply(dims, function(x) any(x == 0L))) &
        !grepl("~", names.gdsn)
        ## what's the pattern with genotype/~data, phase/~data,
        ## annotation/format/DP/~data
    ]
}

.read_gdsdata_sampleInCol <- function(file, node)
{
    ff <- .get_gdsdata_fileFormat(file)
    dims <- .get_gdsdata_dim(file, node)
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

.get_gdsdata_dim <- function(file, node)
{
    f <- openfn.gds(file)
    on.exit(closefn.gds(f))
    
    dim <- objdesp.gdsn(index.gdsn(f, node))$dim
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

.get_gdsdata_dimnames <- function(file, node)
{
    ff <- .get_gdsdata_fileFormat(file)
    dims <- .get_gdsdata_dim(file, node)

    if (ff == "SNP_ARRAY") {
        f <- openfn.gds(file)
        on.exit(closefn.gds(f))
        
        sample.id <- read.gdsn(index.gdsn(f, "sample.id"))
        snp.id <- read.gdsn(index.gdsn(f, "snp.id"))
        rd <- names(get.attr.gdsn(index.gdsn(f, node))) ## ?
        if ("snp.order" %in% rd) {
            dimnames <-
                list(snp.id = as.character(snp.id), sample.id = sample.id)
        } else {
            dimnames <-
                list(sample.id = sample.id, snp.id = as.character(snp.id))
        }
        dimSumm <- lengths(dimnames)
    } else if (ff == "SEQ_ARRAY") {
        f <- seqOpen(file)
        on.exit(seqClose(f))
        
        sample.id <- read.gdsn(index.gdsn(f, "sample.id"))
        variant.id <- read.gdsn(index.gdsn(f, "variant.id"))
        seqSumm <- seqSummary(f, verbose=FALSE)
        dimSumm <- c(
            ploidy = seqSumm$ploidy,
            sample = seqSumm$num.sample,
            variant = seqSumm$num.variant)
        stopifnot(length(variant.id) == dimSumm["variant"])
        stopifnot(length(sample.id) == dimSumm["sample"])
        dimnames <- list(
            ploidy.id = seq_len(dimSumm[1]),
            sample.id = sample.id,
            variant.id = as.character(variant.id)
        )
    }
        ind <- match(dims, dimSumm)
        dimnames <- dimnames[ind]
    if (!is.list(dimnames)) {
        stop(wmsg(
            "The dimnames of GDS dataset '", file,
            "' should be a list!"))
    }
    dimnames
}
## FIXME: if dimnames is shorter than the corresponding dimensions,
## use NULL to extend.??

.read_gdsdata_first_val <- function(file, node){
    dims <- .get_gdsdata_dim(file, node)
    f <- openfn.gds(file)
    on.exit(closefn.gds(f))
    first_val <- readex.gdsn(
        index.gdsn(f, node), sel=as.list(rep(1, length(dims))))
    first_val
}
