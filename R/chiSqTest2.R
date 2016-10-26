#' #' @title Chi square test
#' #' @description Chi square test for IBD sub population data
#' #'
#' #' @param locusI chr and position locus I
#' #' @param locusII chr and position locus II
#' #' @param locusIII chr and position locus III
#' #' @param x ibdData object
#' #' @param spop sub-population index
#' #'
#' #' @importFrom stats pchisq
#' #' @export
#' chiSqTest3Loci.ibd <- function(x, locusI, locusII, locusIII, spop = rep(1, nrow(x))) {
#'   stopifnot(is(x, "ibdData"))
#'   spop <- as.character(spop)
#'   #   label.inx <- sort(as.integer(unique(as.vector(x@ibdData))))
#'   #   npar <- length(unique(unlist(strsplit(x@alleleLabels[label.inx], ""))))
#'   l.list <- lapply(unique(spop), function(i, x) which(i == x), x = spop)
#'   pos <- x@position
#'   chr.all <- unique(pos[1, ])
#'   if (length(chr.all) == 1L)
#'     stop("data are only from one chromosome")
#'   if (locusI[1] == locusII[1])
#'     stop("Loci are at the same chromosome")
#'   inxI <- which(pos[1, ] == locusI[1])
#'   i <- inxI[findInterval(locusI[2], pos[2, inxI])]
#'   inxII <- which(pos[1, ] == locusII[1])
#'   j <- inxII[findInterval(locusII[2], pos[2, inxII])]
#'   inxIII <- which(pos[1, ] == locusIII[1])
#'   k <- inxIII[findInterval(locusIII[2], pos[2, inxIII])]
#'   obs <- 0
#'   exp <- 0
#'   for (l in 1:length(l.list)) {
#'     tables <- tables3V(x@ibdData[l.list[[l]], i],
#'                        x@ibdData[l.list[[l]], j],
#'                        x@ibdData[l.list[[l]], k],
#'                        x@alleleLabels)
#'     obs <- obs + tables[[1]]
#'     exp <- exp + tables[[2]]
#'   }
#'   ndim <- as.integer(round(length(obs)^(1/3), 0))
#'   fullcellinx <- which(exp != 0)
#'   chi <- sum((obs[fullcellinx] - exp[fullcellinx])^2 / exp[fullcellinx])
#'   dfs <- (ndim - 1)^3
#'   p.val <- pchisq(chi, dfs, lower.tail = FALSE)
#'   out <- list(observed = array(obs, c(ndim, ndim, ndim)),
#'               expected = array(exp, c(ndim, ndim, ndim)),
#'               pValue = p.val)
#'   out
#' }
