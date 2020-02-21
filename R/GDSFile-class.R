###----------------
### GDSFile class
###----------------
#' @exportClass GDSFile
#' @rdname GDSFile-class
#' @description \code{GDSFile}: \code{GDSFile} is a light-weight class
#'     to represent a GDS file. It has the `$` completion method to
#'     complete any possible gds nodes. If the slot of `current_path`
#'     in `GDSFile` object represent a valid gds node, it will return
#'     the `GDSArray` of that node directly. Otherwise, it will return
#'     the `GDSFile` object with an updated `current_path`.
#' @aliases GDSFile-class

setClass(
    "GDSFile",
    slots=c(
        file = "character",
        current_path = "character"
    )
)

setMethod("show", "GDSFile", function(object) {
    nodes <- gdsnodes(object)
    nodes <- nodes[startsWith(nodes, object@current_path)]
    cat(
        "class: ", class(object), "\n",
        "file: ", object@file, "\n",
        "current node: ", object@current_path, "\n",
        "subnodes:\n  ", paste(nodes, collapse="\n  "), "\n",
        sep = ""
    )
})

###--------------
### constructor
###--------------

#' GDSFile constructor and methods. 
#' 
#' @name GDSFile
#' @description \code{GDSFile}: the \code{GDSFile} class constructor.
#' @param file the GDS file path.
#' @param current_path the current path to the closest gds node.
#' @export
#' @rdname GDSFile-class
#' @aliases GDSFile-constructor
GDSFile <- function(file, current_path="")
{
    new("GDSFile", file = file, current_path = current_path)
}

###------------
### accessors
###------------

#' @exportMethod gdsfile
#' @rdname GDSFile-class
#' @aliases GDSFile-method
#' @description \code{gdsfile}: \code{file} slot getter for
#'     \code{GDSFile} object.
#' @param object \code{GDSFile} object.
#' @return \code{gdsfile}: the file path of corresponding \code{GDSfile} object.

setMethod("gdsfile", "GDSFile", function(object) object@file)

#' @exportMethod "gdsfile<-"
#' @rdname GDSFile-class
#' @aliases GDSFile-method GDSFile,gdsfile-method
#' @description \code{gdsfile<-}: \code{file} slot setter for
#'     \code{GDSFile} object.
#' @param value the new gds file path
#' @examples
#' file <- SeqArray::seqExampleFileName("gds")
#' gf <- GDSFile(file)
#' gdsfile(gf)
#' 
setReplaceMethod("gdsfile", "GDSFile", function(object, value) {
    new_filepath <- tools::file_path_as_absolute(value)
    BiocGenerics:::replaceSlots(object, file=value, check=FALSE)
})

###--------------------
### dollar completion
###--------------------
#' @importFrom utils .DollarNames
#' @export
.DollarNames.GDSFile <- function(x, pattern = "") {
    nodes <- gdsnodes(x)
    nodes <- nodes[startsWith(nodes, x@current_path)]
    completions <- sub(sprintf("^%s/", x@current_path), "", nodes)
    sub("/.*", "", completions)
}

#' @exportMethod $
#' @rdname GDSFile-class
#' @aliases GDSFile-method
#' @param name the name of gds node
#' @return \code{$}: a \code{GDSFile} with updated \code{@current_path}, or
#'     \code{GDSArray} object if the \code{current_path} is a valid
#'     gds node.
#' 
setMethod("$", "GDSFile", function(x, name)
{
    if (nzchar(x@current_path)) {
        name <- paste(x@current_path, name, sep="/")
    }
    ## check if exist
    nodes <- gdsnodes(x)
    nodes <- nodes[startsWith(nodes, name)]
    pattern <- sprintf("(%s*)/.*$", name)

    if (name %in% sub(pattern, "\\1", nodes)) {
        x@current_path <- name
    } else {
        stop(wmsg("the gds path of '", name, "' does not exist"))
    }
    if (x@current_path %in% gdsnodes(x@file))
        GDSArray(gdsfile(x), x@current_path)
    else
        x
})

###------------
### methods
###------------

#' @exportMethod gdsnodes

setGeneric("gdsnodes", function(x) standGeneric(x), signature="x")

#' @name gdsnodes
#' @rdname GDSFile-class
#' @aliases GDSFile-method gdsnodes,ANY-method gdsnodes,GDSFile-method
#' @description \code{gdsnodes}: to get the available gds nodes from
#'     the \code{GDSFile} object or the file path with extension of
#'     ".gds".
#' @param x a \code{GDSFile} object. or GDS file path (for
#'     \code{gdsnodes()}).
#' @return \code{gdsnodes}: a character vector for the available gds
#'     nodes. When input is GDS file path, it returns all available
#'     gds nodes within the GDS file, no matter there is value or
#'     not. When input is \code{GDSFile} object, it returns only the
#'     gds nodes that could construct unique \code{GDSArray} objects,
#'     which means that the gds node has non-zero-dimensions, and is
#'     actually array, and all \code{GDSArray}s returned from these
#'     nodes are unique (by excluding the gds nodes that has `code{~}
#'     prefix).
#' @examples
#' file <- SNPRelate::snpgdsExampleFileName()
#' gdsnodes(file)
#' file1 <- SeqArray::seqExampleFileName("gds")
#' gdsnodes(file1)
#' gf <- GDSFile(file)
#' gdsnodes(gf)
#' gdsfile(gf)
setMethod("gdsnodes", "ANY", function(x)
{
    ff <- .get_gds_fileFormat(x)
    if (ff == "SEQ_ARRAY") {
        cmnnodes <- c("sample.id", "variant.id", "position",
                         "chromosome", "allele", "genotype",
                         "annotation/id", "annotation/qual",
                         "annotation/filter")
        f <- seqOpen(x)
        on.exit(seqClose(f))
        info.nodes <- ls.gdsn(index.gdsn(f, "annotation/info"))
        format.nodes <- ls.gdsn(index.gdsn(f, "annotation/format"))
        names.gdsn <- c(cmnnodes, paste0("annotation/info/", info.nodes),
                        paste0("annotation/format/", format.nodes))
    } else {
        f <- openfn.gds(x)
        on.exit(closefn.gds(f))
        names.gdsn <- ls.gdsn(f)
        repeat {
            a <- lapply(names.gdsn, function(x) ls.gdsn(index.gdsn(f, x)))
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
    }
    names.gdsn
})

#' @exportMethod gdsnodes
setMethod("gdsnodes", "GDSFile", function(x)
{
    nodes <- gdsnodes(gdsfile(x))
})

## #' @exportMethod gdsnodes
## setMethod("gdsnodes", "GDSFile", function(x)
## {
##     nodes <- gdsnodes(gdsfile(x))
##     isarray <- vapply(nodes, function(node)
##         .get_gdsnode_isarray(gdsfile(x), node),
##         logical(1))
##     permuted_same_data <- grepl("~", nodes)
##     dims <- lapply(nodes, function(node)
##                    .get_gdsnode_dim(gdsfile(x), node))
##     isnull_dims <- vapply(dims, is.null, logical(1))
##     any_zero_dims <- vapply(dims, function(dim) any(dim == 0), logical(1))

##     nodes[isarray & !permuted_same_data & !isnull_dims & !any_zero_dims]
## })


## todo: 1. GDSFile(seqfile) debug.  -- done! 
## 2. remove ".get_gdsnode_isarray". Remove, Check, Documentation.
## 3. remove ".get_gdsnode_first_val"?? used in .extract_array now. 
## 4. check DDF, VE.
