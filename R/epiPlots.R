


#' @importFrom graphics plot lines
#' @export
chr.plot <- function (chr.sizes, chr.radius=.8, gap.rate=.01, add=FALSE, ...) {
  chr.sr <- chr.pos(chr.sizes, gap.rate)
  s.s.chr <- chr.sr[[1]]
  a <- chr.sr[[2]]
  # Edit Chr positio
  if (!add)
    plot(-1:1, -1:1, type="n", xaxs="i", yaxs="i", 
         axes=FALSE, xlab="", ylab="", asp=1)
  for (i in 1:length(chr.sizes)) {
    chr.x <- seq(s.s.chr[i, 1], s.s.chr[i, 2], length.out=100)
    chr.y <- rep(chr.radius, 100)
    chr.x.circ <- chr.y*sin(a*(chr.x-1))
    chr.y.circ <- chr.y*cos(a*(chr.x-1))
    lines(x=chr.x.circ, chr.y.circ, ...)
  }
}



#' @importFrom graphics lines
#' @export
add.lines <- function (chr, pos, chr.sizes, chr.radius=.8, 
                       gap.rate=.01, ...) {
  chr.sr <- chr.pos(chr.sizes, gap.rate)
  s.s.chr <- chr.sr[[1]]
  a <- chr.sr[[2]]
  # Edit aditive effects
  pos <- seq(min(pos), max(pos), length.out=20)
  x  <- pos + s.s.chr[chr, 1]
  x.circ <- chr.radius*sin(a*(x-1))
  y.circ <- chr.radius*cos(a*(x-1))
  lines(x.circ, y.circ, ...)
}



#' @importFrom graphics polygon
#' @export
epi.lines <- function (chr1, pos1, chr2, pos2, chr.sizes, chr.radius=.8, 
                       gap.rate=.01, ...) {
  chr.sr <- chr.pos(chr.sizes, gap.rate)
  s.s.chr <- chr.sr[[1]]
  a <- chr.sr[[2]]
  # Edit epistatic effects
  pos1 <- seq(min(pos1), max(pos1), length.out=20)
  pos2 <- seq(max(pos2), min(pos2), length.out=20)
  x1  <- pos1+s.s.chr[chr1, 1]
  x2  <- pos2+s.s.chr[chr2, 1]
  x1.circ <- chr.radius*sin(a*(x1-1))
  y1.circ <- chr.radius*cos(a*(x1-1))
  x2.circ <- chr.radius*sin(a*(x2-1))
  y2.circ <- chr.radius*cos(a*(x2-1))
  x1x2 <- c(x1.circ[1], 0, x2.circ[1])
  y1y2 <- c(y1.circ[1], 0, y2.circ[1])
  x2x1 <- c(x2.circ[length(x2.circ)], 0, x1.circ[length(x1.circ)]) 
  y2y1 <- c(y2.circ[length(y2.circ)], 0, y1.circ[length(y1.circ)]) 
  x1x2.poly <- y1y2.poly <- x2x1.poly <- y2y1.poly <- double(100)
  x1x2.poly[c(1, 100)] <- x1x2[c(1, 3)]
  y1y2.poly[c(1, 100)] <- y1y2[c(1, 3)]
  x2x1.poly[c(1, 100)] <- x2x1[c(1, 3)]
  y2y1.poly[c(1, 100)] <- y2y1[c(1, 3)]
  for (i in 2:99) {
    x1x2i.poly <- y1y2i.poly <- x2x1i.poly <- y2y1i.poly <- 0
    const <- (1-i*.01)^2
    for (j in 1:3) {
      x1x2i.poly <- x1x2i.poly+const*x1x2[j]
      y1y2i.poly <- y1y2i.poly+const*y1y2[j]
      x2x1i.poly <- x2x1i.poly+const*x2x1[j]
      y2y1i.poly <- y2y1i.poly+const*y2y1[j]
      const <- const*(3-j)/j*i*.01/(1-i*.01)
    }
    x1x2.poly[i] <- x1x2i.poly
    y1y2.poly[i] <- y1y2i.poly
    x2x1.poly[i] <- x2x1i.poly
    y2y1.poly[i] <- y2y1i.poly
  }
  polygon(c(x1x2.poly, x2.circ, x2x1.poly, x1.circ), 
          c(y1y2.poly, y2.circ, y2y1.poly, y1.circ), ...) 
}



chr.pos <- function(chr.sizes, gap.rate) {
  gap.chr <- round(sum(chr.sizes)*gap.rate, 0)
  i.chr <- 0
  s.s.chr <- matrix(NA, length(chr.sizes), 2)
  for (i in 1:length(chr.sizes)) {
    i.chr <- i.chr+gap.chr 
    s.s.chr[i, 1] <- i.chr+1
    i.chr <- i.chr+chr.sizes[i]
    s.s.chr[i, 2] <- i.chr
    i.chr <- i.chr+gap.chr 
  }
  a <- pi*2/i.chr
  out <- list(s.s.chr, a)
  out
}
