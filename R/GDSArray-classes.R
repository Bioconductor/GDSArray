#' GDSArraySeed / GDSArray classes and constructors for gds data.

### -----------------------------------------------------
### GDSArraySeed class 
###

#' @importClassesFrom DelayedArray Array
#' @export
#' @rdname GDSArray-classes
setClass("GDSArraySeed",
         contains = "Array", ## from DelayedArray: A virtual class with no slots
                             ## to be extended by concrete subclasses with
                             ## an array-like semantic.
    slots = c(
        file="character",   # Absolute path to the gds file so the object
                            # doesn't break when the user changes the working
                            # directory (e.g. with setwd()).
        name="character",   # Name of the dataset in the gds file.
        dim = "integer",
        dimnames = "list",
        permute = "logical",
        first_val = "ANY"      ## remove this slot. 
    )
)

###
## show method for GDSArraySeed object
###
setMethod(
    "show", "GDSArraySeed",
    function(object) {
        cat("GDSArraySeed\n",
            "gds file: ", object@file, "\n",
            "array data: ", object@name, "\n",
            ## "dim: ", nrow(object), " x ", ncol(object), "\n",
            "dim: ", paste(dim(object), collapse=" x "), "\n",
            sep="")
    }
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### extract_array()
###
#' @importMethodsFrom DelayedArray extract_array
.extract_array_from_GDSArraySeed <- function(x, index)
{
    ans_dim <- DelayedArray:::get_Nindex_lengths(index, dim(x))
    if (any(ans_dim == 0L)){
        ans <- x@first_val[0]  ## integer(0) / character(0)
        dim(ans) <- ans_dim
    } else {
        f <- openfn.gds(x@file)
        on.exit(closefn.gds(f))
        if(x@permute){
            permdim <- rev(seq_len(length(index)))
            index <- index[permdim] ## multi-dimensional supported
            ans <- aperm(readex.gdsn(index.gdsn(f, x@name), index), permdim)
            ## same: ans <- aperm(readex.gdsn(index.gdsn(f, x@name), index))
        } else {
            ans <- readex.gdsn(index.gdsn(f, x@name), index)
        }
    }
    ans
}
setMethod("extract_array", "GDSArraySeed", .extract_array_from_GDSArraySeed)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### GDSArraySeed constructor
###

#' GDSArraySeed
#' @description \code{GDSArraySeed}: The function to generate a
#'     GDSArraySeed for the later conversion from gds file into
#'     GDSArray.
#' @param file the gds file name.
#' @param name the gds array node to be read into GDSArraySeed / GDSArray. For
#'     \code{GDSArray}, the default value for \code{name} is the
#'     genotype data.
#' @export
#' @rdname GDSArray-classes
#' @importFrom SNPRelate snpgdsOpen snpgdsClose
#' @importFrom SeqArray seqOpen seqClose seqSummary
#' @import methods
#' @import gdsfmt
#' @importFrom tools file_path_as_absolute

GDSArraySeed <- function(file, name=NA)
{
    if (!isSingleString(file))
        stop(wmsg("'file' must be a single string specifying the path to ",
                  "the gds file where the dataset is located."))
    if (!isSingleStringOrNA(name))
        stop("'type' must be a single string or NA")
    file <- file_path_as_absolute(file)

    ## ff <- .get_gdsdata_fileFormat(file)

    ## if (ff == "SNP_ARRAY") {
    ##     f <- snpgdsOpen(file)
    ##     on.exit(snpgdsClose(f))
    ## } else if (ff == "SEQ_ARRAY") {
    ##     f <- seqOpen(file)
    ##     on.exit(seqClose(f))
    ## } else {
    ##     f <- openfn.gds(file)
    ##     on.exit(closefn.gds(f))
    ## }
    
    ## ## check if the node is array data. 
    ## arrayNodes <- .get_gdsdata_arrayNodes(f)
    ## if(!name %in% arrayNodes){
    ##     stop(wmsg("the `name` node must be an array data."))
    ## }

    dims <- .get_gdsdata_dim(file, node = name)
    dimnames <- .get_gdsdata_dimnames(file, node = name)

    if (!identical(lengths(dimnames, use.names=FALSE), dims)) {
        stop(wmsg("the lengths of dimnames is not consistent with data dimensions."))
    }

    first_val <- .read_gdsdata_first_val(file, node = name)

    if (length(dims) == 1)
        permute = FALSE
    else 
        permute = !.read_gdsdata_sampleInCol(file, node = name)

    if (permute) {
        dims <- rev(dims)
        dimnames <- dimnames[rev(seq_len(length(dimnames)))]
    }
    new2("GDSArraySeed", file=file,
         name=name,
         dim=dims,
         dimnames = dimnames,
         permute = permute,
         first_val = first_val
         )
}


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### GDSArray and GDSMatrix objects
###
### We define these classes only for cosmetic reasons i.e. to hide the
### DelayedArray and DelayedMatrix classes from the user. The user will see
### and manipulate GDSArray and GDSMatrix objects instead of DelayedArray
### and DelayedMatrix objects.
###

#' @importClassesFrom DelayedArray DelayedArray DelayedMatrix
setClass("GDSArray", contains="DelayedArray")
setClass("GDSMatrix", contains=c("DelayedMatrix", "GDSArray"))

### Automatic coercion method from GDSArray to GDSMatrix (muted for
### higher dimensions) this function works only when GDSArray is
### 2-dimensional, otherwise, it fails.
setAs("GDSArray", "GDSMatrix", function(from) new("GDSMatrix", from))    
setAs("GDSMatrix", "GDSArray", function(from) from)

### For internal use only.
setMethod("matrixClass", "GDSArray", function(x) "GDSMatrix")

.validate_GDSArray <- function(x)
{
    if (!is(x@seed, "GDSArraySeed"))
        return(wmsg("'x@seed' must be a GDSArraySeed object"))
    if (!DelayedArray:::is_pristine(x))
        return(wmsg("'x' carries delayed operations"))
    TRUE
}

#' @importFrom S4Vectors setValidity2
setValidity2("GDSArray", .validate_GDSArray)

## setAs("ANY", "GDSMatrix",
##       function(from) as(as(from, "GDSArray"), "GDSMatrix")
##     )


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor
###

#' @importFrom DelayedArray DelayedArray
setMethod("DelayedArray", "GDSArraySeed",
          function(seed) DelayedArray:::new_DelayedArray(seed, Class="GDSArray")
          )

#' @description \code{GDSArray}: The function to convert a gds file
#'     into the GDSArray data structure.
#' @param index the class of the index slot of the
#'     GDSArray. \code{IndexList} is a reference class list defined in
#'     package \code{DelayedArray}, in order to have different
#'     GDSArrays to share the index slot, for a uniform subsetting
#'     performance.
#' @export
#' @rdname GDSArray-classes
#' @examples
#' file <- SNPRelate::snpgdsExampleFileName()
#' allnodes <- gdsNodes(file)  ## print all available gds nodes in file.
#' allnodes
#' GDSArraySeed(file, "genotype")
#' GDSArray(file)
#' GDSArray(GDSArraySeed(file, "genotype"))
#'
#' file1 <- SeqArray::seqExampleFileName("gds")
#' allnodes1 <- gdsNodes(file1)  ## print all available gds nodes in file1. 
#' allnodes1
#' GDSArraySeed(file1, "genotype/data")
#' GDSArray(file1)
#' GDSArray(file1, "variant.id")
#' GDSArray(file1, "sample.annotation/family")
#' GDSArray(file1, "annotation/format/DP/data")
#' GDSArray(file1, "annotation/info/DP")

GDSArray <- function(file, name=NA, index=c("list", "IndexList")){
    if (is(file, "GDSArraySeed")) {
        seed <- file
    } else {
        ff <- .get_gdsdata_fileFormat(file)
        if (ff == "SNP_ARRAY") {
            if (is.na(name)) name <- "genotype"
        } else if (ff == "SEQ_ARRAY") {
            if (is.na(name)) name <- "genotype/data"
        }
        seed <- GDSArraySeed(file, name)
    }
    to <- "GDSArray"
    if (length(dim(seed)) == 2L)
        to <- "GDSMatrix"
    as(DelayedArray(seed), to)
}

