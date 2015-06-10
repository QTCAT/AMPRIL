
#' @title Chi square test
#' @description Chi square test for IBD sub population data
#' @param x ibdData object 
#' @param spop sub-population index
#' @export
chiSqTest.ibd <- function (x, spop = rep(1, nrow(x))) {
  stopifnot(is(x, "ibdData"))
  spop <- as.character(spop)
#   label.inx <- sort(as.integer(unique(as.vector(x@ibdData))))
#   npar <- length(unique(unlist(strsplit(x@alleleLabels[label.inx], ""))))
  l.list <- lapply(unique(spop), function(i, x) which(i == x), x=spop)
  chr.all <- unique(x@position[1, ])
  if (length(chr.all) == 1L) 
    stop("data are only from one chromosome")
  n.chrMin <- sum(min(chr.all) == x@position[1, ])
  n.chrMax <- sum(max(chr.all) == x@position[1, ])
  out <- matrix(NA, ncol(x)-n.chrMax, ncol(x)-n.chrMin)
  for (chr in chr.all[-length(chr.all)]) {
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
        dfs <- (length(nonzc)-1) * (length(nonzr)-1)
        out.j[[j]] <- pchisq(chi, dfs, lower.tail=FALSE)
      }
      out.i[[i]] <- unlist(out.j)
    }
    out[i.inx, j.inx-n.chrMin] <- do.call("rbind", out.i)
  }
  rownames(out) <- colnames(x)[1:nrow(out)]
  colnames(out) <- colnames(x)[-1:-n.chrMin]
  out
} # chiSqTest.ibd
