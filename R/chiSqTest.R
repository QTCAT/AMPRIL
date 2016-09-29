#' @title Chi square test
#' @description Chi square test for IBD sub population data
#' 
#' @param x ibdData object 
#' @param spop sub-population index
#' 
#' @importFrom methods is
#' @export
chiSqTest.ibd <- function(x, spop = rep(1, nrow(x))) {
  stopifnot(is(x, "ibdData"))
  spop <- as.character(spop)
  l.list <- lapply(unique(spop), function(i, x) which(i == x), x = spop)
  chr.all <- unique(x@position[1, ])
  if (length(chr.all) == 1L) 
    stop("data are only from one chromosome")
  out <- list()
  for (chr in chr.all[-length(chr.all)]) { # chr <- 1
    i.inx <- which(chr == x@position[1, ])
    j.inx <- which(chr < x@position[1, ])
    out.i <- list()
    for (i in i.inx) { 
      out.j <- list()
      for (j in j.inx) { 
        obs <- 0
        exp <- 0
        for (l in 1:length(l.list)) { 
          tables <- chiSqTables(x@ibdData[l.list[[l]], i], 
                                x@ibdData[l.list[[l]], j], 
                                x@alleleLabels)
          obs <- obs + tables[[1]]
          exp <- exp + tables[[2]]
        }
        obs <- matrix(obs, sqrt(length(obs)))
        exp <- matrix(exp, sqrt(length(exp)))
        nonzc <- which(colSums(exp) > 0)
        nonzr <- which(rowSums(exp) > 0)
        obs <- obs[nonzr, nonzc]
        exp <- exp[nonzr, nonzc]
        chi <- sum((obs - exp)^2 / exp)
        dfs <- (length(nonzc) - 1) * (length(nonzr) - 1)
        out.j[[j]] <- pchisq(chi, dfs, lower.tail = FALSE)
      }
      out.i[[i]] <- unlist(out.j)
    }
    out.ij <- data.frame(chr1 = rep(x@position[1, i.inx], each = length(j.inx)),
                         pos1 = rep(x@position[2, i.inx], each = length(j.inx)),
                         chr2 = rep(x@position[1, j.inx], times = length(i.inx)),
                         pos2 = rep(x@position[2, j.inx], times = length(i.inx)),
                         pValues = unlist(out.i))
    rownames(out.ij) <- paste(rep(colnames(x)[i.inx], each = length(j.inx)), 
                 rep(colnames(x)[j.inx], times = length(i.inx)), 
                 sep = ":")
    out[[chr]] <- out.ij
  }
  out <- do.call("rbind", out)
  class(out) <- c("epiScan", class(out))
  out
}


#' @title Chi square test
#' @description Chi square test for IBD sub population data
#' 
#' @param locusI chr and position locus I
#' @param locusII chr and position locus II
#' @param x ibdData object 
#' @param spop sub-population index
#' 
#' @importFrom stats pchisq
#' @export
chiSqTest2Loci.ibd <- function (x, locusI, locusII, spop = rep(1, nrow(x))) {
  stopifnot(is(x, "ibdData"))
  spop <- as.character(spop)
  #   label.inx <- sort(as.integer(unique(as.vector(x@ibdData))))
  #   npar <- length(unique(unlist(strsplit(x@alleleLabels[label.inx], ""))))
  l.list <- lapply(unique(spop), function(i, x) which(i == x), x = spop)
  pos <- x@position
  chr.all <- unique(pos[1, ])
  if (length(chr.all) == 1L) 
    stop("data are only from one chromosome")
  if (locusI[1] == locusII[1]) 
    stop("Loci are at the same chromosome")
  inxI <- which(pos[1, ] == locusI[1])
  i <- inxI[findInterval(locusI[2], pos[2, inxI])]
  inxII <- which(pos[1, ] == locusII[1])
  j <- inxII[findInterval(locusII[2], pos[2, inxII])]
  obs <- 0
  exp <- 0
  for (l in 1:length(l.list)) { 
    tables <- chiSqTables(x@ibdData[l.list[[l]], i], 
                          x@ibdData[l.list[[l]], j], 
                          x@alleleLabels)
    obs <- obs+tables[[1]]
    exp <- exp+tables[[2]]
  }
  obs <- matrix(obs, sqrt(length(obs)))
  exp <- matrix(exp, sqrt(length(exp)))
  nonzc <- which(colSums(exp) > 0)
  nonzr <- which(rowSums(exp) > 0)
  obs <- obs[nonzr, nonzc]
  exp <- exp[nonzr, nonzc]
  chi <- sum((obs - exp)^2 / exp)
  dfs <- (length(nonzc) - 1) * (length(nonzr)-1)
  p.val <- pchisq(chi, dfs, lower.tail = FALSE)
  out <- list(observed = obs,
              expected = exp,
              p.value = p.val)
  out
} # chiSqTest.ibd
