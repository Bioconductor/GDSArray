
setClass(
    "GDSlight",
    slots=c(
        file = "character",
        current_path = "character"
    )
)

setMethod("show", "GDSlight", function(object) {
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

#' GDSlight constructor and methods. 
#' 
#' @name GDSlight
#' @description GDSlight is a light-weight class to represent a GDS
#'     file. It has the \code{$} method defined as a convenient input
#'     for \code{GDSArray} constructor.
#' @export
#' @rdname GDSlight-class
#' @aliases GDSlight-constructor
GDSlight <- function(file, current_path="")
{
    new("GDSlight", file = file, current_path = current_path)
}

.DollarNames.GDSlight <- function(x, pattern = "") {
    nodes <- gdsNodes(x@file)
    nodes <- nodes[startsWith(nodes, x@current_path)]
    completions <- sub(sprintf("^%s/", x@current_path), "", nodes)
    sub("/.*", "", completions)
}

#' @exportMethod $
#' @rdname GDSlight-class
#' @aliases GDSlight-method
#' @details \code{$} takes \code{GDSlight} object as input, and output
#'     by updating \code{@current_path}.
setMethod("$", "GDSlight", function(x, name)
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
