#' @title read snp data as ibdData object
#' @description read snp data as ibdData object
#' 
#' @param file the name of the file which the data are to be read from.
#' @param header a logical value indicating whether the file contains the 
#' names of the variables as its first line
#' @param sep the field separator character. Values on each line of the file 
#' are separated by this character.
#' @param quote the set of quoting characters. To disable quoting altogether, 
#' use quote = "".
#' @param row.names a number giving the column of the table which contains 
#' the row names.
#' @param position.ncol integer value indicating whether the how many columns 
#' contain position information.
#' @param na.strings a character vector of strings which are to be interpreted 
#' as NA values. 
#' @param nrows integer: the maximum number of rows to read in.
#' 
#' @importFrom methods new
#' @export
read.ibdData <- function(file, header=TRUE, sep=" ", quote='"', 
                         row.names=FALSE, position.ncol=3L, 
                         na.strings="NA", nrows=-1) {
  temp <- read_ibdData(file, header, sep, quote, 
                       row.names, position.ncol, na.strings, nrows)
  if (identical(temp$indivNames, character(0))) {
    temp$indivNames <- paste0("indiv", seq_len(nrow(temp$ibdData)))
  }
  if (identical(temp$lociNames, character(0))) {
    temp$lociNames <- paste0("loci", seq_len(ncol(temp$ibdData)))
  }
  if (identical(temp$position, integer(0))) {
    temp$position <- NULL
  }
  out <- new("ibdData",
             ibdData=temp$ibdData,
             indivNames=temp$indivNames,
             lociNames=temp$lociNames,
             alleleLabels=temp$alleleLabels,
             position=temp$position)
  out
} # read.ibdData


#' @title Sub ibdData
#' @description ...
#' 
#' @param x ibdData object
#' @param i indices specifying elements to extract or replace. Indices are n
#' umeric or character vectors
#' @param j indices specifying elements to extract or replace. Indices are 
#' numeric or character vectors
#' @param ... not implemented
#' @param drop not implemented
#' 
#' @importFrom methods setMethod signature new
#' @export
setMethod("[", signature(x="ibdData", i="ANY", j="ANY", drop="missing"),
          function(x, i, j, ..., drop) {
            if (!missing(i)) {
              if (is.character(i)) {
                i <- match(i, x@indivNames)
              }
            }
            if (!missing(j)) {
              if (is.character(j)) {
                j <- match(j, x@lociNames)
              }
            }
            out <- new("ibdData",
                       ibdData=x@ibdData[i, j, drop=FALSE],
                       indivNames=x@indivNames[i],
                       lociNames=x@lociNames[j],
                       alleleLabels=x@alleleLabels,
                       position=if(is.null(x@position)) NULL 
                                else x@position[, j, drop=FALSE])
            out
          }
) # `[`

#' @title ibdData as matrix
#' @description ...
#' 
#' @param x ibdData object
#' @param ... not implemented
#' 
#' @importFrom methods setMethod signature
#' @export
setMethod("as.matrix", signature(x="ibdData"),
          function(x, ...) {
            at <- attributes(x@ibdData)
            out <- c(NA, x@alleleLabels)[as.integer(x@ibdData)+1L]
            attributes(out) <- at
            colnames(out) <- x@lociNames
            rownames(out) <- x@indivNames
            attr(out, "position") <- x@position
            out
          }
) # as.matrix


#' @title Get dimnames
#' @description Retrieve or set the dimnames of an object.
#' 
#' @param x ibdData object
#' 
#' @importFrom methods setMethod signature
#' @export
setMethod("dimnames", signature(x="ibdData"),
          function (x) {
            out <- list(x@indivNames, x@lociNames)
            out
          }
) # dimnames


# #'  @title Assign dimnames
# #'  @description Retrieve or set the dimnames of an object.
# #'
# #'  @param x ibdData object
# #'  @param value  List of names
# #'  @details ...
# #'
# #'  @importFrom methods setMethod signature
# #' @export
# setMethod("dimnames<-", signature(x="ibdData", value="list"),
#           function(x, value) {
#             d <- dim(x)
#             if (!is.list(value) || length(value) != 2L ||
#                   !(is.null(v1 <- value[[1L]]) || length(v1) == d[1L]) ||
#                   !(is.null(v2 <- value[[2L]]) || length(v2) == d[2L]))
#               stop(gettextf("invalid dimnames given for %s object",
#                             dQuote(class(x))), domain=NA)
#             if(!is.null(v1)) x@indivNames <- as.character(v1)
#             if(!is.null(v2)) x@lociNames <- as.character(v2)
#             x
#           }
# ) # dimnames


#' @title Get dim
#' @description ..
#' 
#' @param x ibdData object
#' 
#' @importFrom methods setMethod signature
#' @export
setMethod("dim", signature(x="ibdData"),
          function(x) {
            out <- dim(x@ibdData)
            out
          }
) # dim


#' @title Frequency of NAs in ibdData
#' @description ...
#' 
#' @param object ibdData object
#' 
#' @export
na.freq <- function (object) {
  temp <- apply(object@ibdData, 2, function(x){sum(x == 0x00)})
  out <- temp/nrow(object@ibdData)
  names(out) <- object@lociNames
  out
} # "na.freq"


#' @title Get position from ibdData
#' @description ...
#' 
#' @param object ibdData object
#' 
#' @export
getPos <- function (object) {
  out <- object@position
  if (!is.null(out)) {
    colnames(out) <- object@lociNames
  } else {
    cat("No position information available")
  }
  out
} # getPos


#' @title Allele frequency
#' @description ...
#' 
#' @param x ibdData object
#' 
#' @export
allele.freq <- function (x) {
  alleles <- unique(unlist(strsplit(x@alleleLabels, "")))
  out <- apply(x@ibdData, 2, 
               function(x, aa, a){
                 xa <- unlist(strsplit(aa[as.integer(x)], ""))
                 xa <- factor(xa, a)
                 table(xa)
               }, aa=x@alleleLabels, a=alleles)
  out <- out/(2*nrow(x))
  colnames(out) <- colnames(x)
  out
} # allele.freq


#' @title Heterozygosity
#' @description ...
#' 
#' @param x ibdData object
#' @param dim interger for dimention
#' 
#' @export
het.freq <- function (x, dim=c(1, 2)) {
  dim <- dim[1]
  if (nrow(x@position) == 3 & dim == 1L) {
    feq.scale <- x@position[3, ]-x@position[2, ]
  } else {
    feq.scale <- rep(1L, dim(x)[dim])
  }
  labelFreq <- sapply(strsplit(x@alleleLabels, ""), 
                      function(x) length(unique(x)))-1L
  hetFreq <- apply(x@ibdData, dim, 
                   function(x, l, s) sum(l[as.integer(x)]*s)/sum(s), 
                   l=labelFreq, s=feq.scale)
  if (dim == 1) {
    names(hetFreq) <- x@indivNames
  }
  if (dim == 2) {
    names(hetFreq) <- x@lociNames
  }
  hetFreq
} # het.freq
