
setClass(
    "GDSFile",
    slots=c(
        file = "character",
        current_path = "character"
    )
)

setMethod("show", "GDSFile", function(object) {
    nodes <- gdsNodes(object@file)
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
#' @description GDSFile is a light-weight class to represent a GDS
#'     file. It has the \code{$} method defined as a convenient input
#'     for \code{GDSArray} constructor.
#' @export
#' @rdname GDSFile-class
#' @aliases GDSFile-constructor
GDSFile <- function(file, current_path="")
{
    new("GDSFile", file = file, current_path = current_path)
}

.DollarNames.GDSFile <- function(x, pattern = "") {
    nodes <- gdsNodes(x@file)
    nodes <- nodes[startsWith(nodes, x@current_path)]
    completions <- sub(sprintf("^%s/", x@current_path), "", nodes)
    sub("/.*", "", completions)
}

#' @exportMethod $
#' @rdname GDSFile-class
#' @aliases GDSFile-method
#' @details \code{$} takes \code{GDSFile} object as input, and output
#'     by updating \code{@current_path}.
setMethod("$", "GDSFile", function(x, name)
{
    if (nzchar(x@current_path)) 
        name <- paste(x@current_path, name, sep="/")
    x@current_path <- name
    ## check if exist
    if (x@current_path %in% gdsNodes(x@file))
        GDSArray(x)
    else
        x
})

###------------
### accessors
###------------

#' @exportMethod gdsfile
#' @rdname GDSFile-class
#' @aliases GDSFile-method
#' @description \code{file} slot getter for \code{GDSFile} object.
setMethod("gdsfile", "GDSFile", function(x) x@file)

###------------
### methods
###------------

#' @exportMethod gdsNodes
#' @name gdsNodes
#' @rdname GDSFile-class
#' @aliases GDSFile-method
#' @description to get the available gds nodes for the \code{GDSFile} object.
#' @return a character vector for the available gds nodes.
setGeneric("gdsNodes", function(x) standGeneric(x))
setMethod("gdsNodes", "GDSFile", function(x) gdsNodes(gdsfile(GDSFile)) )
