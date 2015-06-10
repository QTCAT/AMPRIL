
#' @title Class "ibdData" for a IBD-block matrix (S4 class)
#' @name ibdData-class
#' @docType class
#' @description A S4 class for a SNP matrix. Storing SNP information using a 
#' byte-level (raw) storage scheme jointly with positions.
#' @slot ibdData A matrix of IBD-block stored in type 'raw'.
#' @slot indivNames A vector of individual names.
#' @slot lociNames A vector of loci names.
#' @slot alleleLabels A vector of the block-allele label names.
#' @slot position A matrix with three rows. The first row contains the 
#' chromosomes, the second row contains the start positions of a block, and the 
#' third row contains the end population of a block.
#' @importFrom methods setClass
#' @export
setClass(Class="ibdData",
         slots=c(ibdData="matrix",
                 indivNames="character", 
                 lociNames="character",
                 alleleLabels="character",
                 position="matrix"))
