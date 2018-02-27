.get_gdsdata_fileFormat <- function(file)
{
    f <- openfn.gds(file)
    on.exit(closefn.gds(f))
    ff <- get.attr.gdsn(f$root)$FileFormat
    ff
}

.get_gdsdata_allNodes <- function(gdsfile)
{
    stopifnot(inherits(gdsfile, "gds.class"))
    names.gdsn <- ls.gdsn(gdsfile)
    repeat {
        a <- lapply(names.gdsn, function(x) ls.gdsn(index.gdsn(gdsfile, x)))
        if (all(lengths(a)==0L)) {
            break
        } else {
            a[lengths(a)==0] <- ""
        n <- rep(names.gdsn, lengths(a))
        all.gdsn <- paste(n, unlist(a), sep="/")
        all.gdsn <- sub("/$", "", all.gdsn)
        names.gdsn <- all.gdsn
        }
    }
    names.gdsn
}

## array data with >1 dimensions to pass into assays(se)
.get_gdsdata_arrayNodes <- function(gdsfile)
{
    stopifnot(inherits(gdsfile, "gds.class"))
    names.gdsn <- .get_gdsdata_allNodes(gdsfile)
    isarray <- vapply(names.gdsn, function(x) {
        objdesp.gdsn(index.gdsn(gdsfile, x))$is.array
    }, logical(1))
    dims <- lapply(names.gdsn, function(x)objdesp.gdsn(index.gdsn(gdsfile, x))$dim)
    ## names(dims) <- all.gdsn
    names.gdsn[isarray &
             lengths(dims) > 1 &
             !unlist(lapply(dims, function(x) any(x == 0L))) &
             !grepl("~", names.gdsn)
             ## what's the pattern with genotype/~data, phase/~data,
             ## annotation/format/DP/~data
             ]
}

.read_gdsdata_sampleInCol <- function(gdsfile, node, fileFormat)
{
    stopifnot(inherits(gdsfile, "gds.class"))
    if (fileFormat == "SNP_ARRAY") {
        rd <- names(get.attr.gdsn(index.gdsn(gdsfile, node)))
        if ("snp.order" %in% rd) sampleInCol <- TRUE   ## snpfirstdim (in row)
        if ("sample.order" %in% rd) sampleInCol <- FALSE
    } else if (fileFormat == "SEQ_ARRAY") {
        seqSumm <- seqSummary(gdsfile, verbose=FALSE)
        dimSumm <- c(ploidy = seqSumm$ploidy,
                     sample = seqSumm$num.sample,
                     variant = seqSumm$num.variant)
        dims <- .get_gdsdata_dim(gdsfile, node)
        ind <- match(dimSumm[c("variant", "sample")], dims)
        if (ind[1] < ind[2]) {
            sampleInCol <- TRUE   ## rewrite for general cases: format/DP.
        } else {
            sampleInCol <- FALSE
        }
    }
    sampleInCol
}

.get_gdsdata_dim <- function(gdsfile, node)
{
    stopifnot(inherits(gdsfile, "gds.class"))
    dim <- objdesp.gdsn(index.gdsn(gdsfile, node))$dim
    if (!is.integer(dim)) {
        if (any(dim > .Machine$integer.max)) {
            dim_in1string <- paste0(dim, collapse=" x ")
            stop(wmsg("The dimensions of GDS dataset '", file, "' are: ",
                      dim_in1string, "\n\nThe GDSArray package only ",
                      "supports datasets with all dimensions <= 2^31-1",
                      " (this is ", .Machine$integer.max, ") at the moment."))
        }
    }
    dim <- as.integer(dim)
    dim
}

.get_gdsdata_dimnames <- function(gdsfile, node, fileFormat)
{
    stopifnot(inherits(gdsfile, "gds.class"))
    dims <- .get_gdsdata_dim(gdsfile, node)
    sample.id <- read.gdsn(index.gdsn(gdsfile, "sample.id"))
    if (fileFormat == "SNP_ARRAY") {
        snp.id <- read.gdsn(index.gdsn(gdsfile, "snp.id"))
        rd <- names(get.attr.gdsn(index.gdsn(gdsfile, node))) ## ?
        if ("snp.order" %in% rd) {
            dimnames <- list(snp.id = as.character(snp.id), sample.id = sample.id)
        } else {
            dimnames <- list(sample.id = sample.id, snp.id = as.character(snp.id))
        }
    } else if (fileFormat == "SEQ_ARRAY") {
        variant.id <- read.gdsn(index.gdsn(gdsfile, "variant.id"))
        seqSumm <- seqSummary(gdsfile, verbose=FALSE)
        dimSumm <- c(ploidy = seqSumm$ploidy,
                     sample = seqSumm$num.sample,
                     variant = seqSumm$num.variant)
        stopifnot(length(variant.id) == dimSumm["variant"])
        stopifnot(length(sample.id) == dimSumm["sample"])
        dimnames <- list(
            ploidy.id = seq_len(dimSumm[1]),
            sample.id = sample.id,
            variant.id = as.character(variant.id)
        )
        ind <- match(dims, dimSumm)
        dimnames <- dimnames[ind]
    }
    if (!is.list(dimnames)) {
        stop(wmsg("The dimnames of GDS dataset '", file, "' should be a list!"))
    }
    dimnames
}
## if dimnames is shorter than the corresponding dimensions, use NULL to extend.??

.read_gdsdata_first_val <- function(gdsfile, node){
    dims <- .get_gdsdata_dim(gdsfile, node)
    first_val <- readex.gdsn(index.gdsn(gdsfile, node), sel=as.list(rep(1, length(dims))))
    first_val
}
