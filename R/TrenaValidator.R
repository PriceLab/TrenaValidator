#----------------------------------------------------------------------------------------------------
#' @import methods
#' @importFrom AnnotationDbi select
#' @import org.Hs.eg.db
#' @import TrenaProject
#' @import TrenaProjectHG38.generic
#' @import MotifDb
#' @import phastCons7way.UCSC.hg38
#' @import GenomicScores
#' @import trena
#'
#' @title TrenaValidator-class
#'
#' @name TrenaValidator-class
#' @rdname TrenaValidator-class
#' @aliases TrenaValidator
#' @exportClass TrenaValidator
#'

.TrenaValidator <- setClass("TrenaValidator",
                            representation=representation(
                               TF="character",
                               targetGene="character",
                               tbl.benchmark="data.frame",
                               state="environment",
                               trenaProject="TrenaProject",
                               quiet="logical"
                               ))

#----------------------------------------------------------------------------------------------------
setGeneric('fallbackEnhancers',  signature='obj', function(obj) standardGeneric('fallbackEnhancers'))
setGeneric('findEnhancers',  signature='obj', function(obj, eliteOnly) standardGeneric('findEnhancers'))
setGeneric('getTFBS',  signature='obj', function(obj, tbl.regions, fimo.threshold, conservation.threshold, pwmFile)
              standardGeneric('getTFBS'))
setGeneric('setMatrix',   signature='obj', function(obj, matrix) standardGeneric('setMatrix'))
setGeneric('buildModel',  signature='obj', function(obj) standardGeneric('buildModel'))
#----------------------------------------------------------------------------------------------------
#' Define an object of class TrenaValidator
#'
#' @description
#' A TF/targetGene/tissue combination
#'
#' @rdname TrenaValidator-class
#'
#' @export
#'
#' @return An object of the TrenaValidator class
#'

TrenaValidator <- function(TF, targetGene, tbl.benchmark, quiet=TRUE)
{
   source("~/github/fimoService/batchMode/fimoBatchTools.R")   # works on hagfish & khaleesi

   state <- new.env(parent=emptyenv())
   state$matrix <- matrix()

   tp <- TrenaProjectHG38.generic()
   setTargetGene(tp, targetGene)

   .TrenaValidator(TF=TF, targetGene=targetGene, tbl.benchmark=tbl.benchmark,
                   state=state, trenaProject=tp, quiet=quiet)

} # TrenaValidator, the constructor
#----------------------------------------------------------------------------------------------------
#' a simple single 10k regulatory region used if no geneHancers values are avaialable
#'
#' @param obj An instance of the TrenaValidator class
#'
#' @return data.frame
#'
#' @export
#'
#' @aliases fallbackEnhancers
#' @rdname fallbackEnhancers

setMethod('fallbackEnhancers',  'TrenaValidator',

   function(obj) {
     tbl.geneInfo <- getTranscriptsTable(obj@trenaProject)
     tbl.fallbackEnhancers <- with(tbl.geneInfo,
                                   data.frame(chrom=chrom, start=tss-5000, end=tss+5000,
                                              stringsAsFactors=FALSE))
     tbl.fallbackEnhancers
     })

#---------------------------------------------------------------------------------------------------
#' query and store the GeneHancer regions for the targetGene
#'
#' @param obj An instance of the TrenaValidator class
#' @param eliteOnly keep all regions, on only those with two attestations
#'
#' @return nothing
#'
#' @export
#'
#' @aliases findEnhancers
#' @rdname findEnhancers

setMethod('findEnhancers',  'TrenaValidator',

      function(obj, eliteOnly){
         tbl.enhancers <- getEnhancers(obj@trenaProject)
         if(eliteOnly)
            tbl.enhancers <- subset(tbl.enhancers, elite)
         if(nrow(tbl.enhancers) == 0)
            return(fallbackEnhancers(obj))
         return(tbl.enhancers)
         })

#----------------------------------------------------------------------------------------------------
#' return the TFBS
#'
#' @param obj An instance of the TrenaValidator class
#'
#' @return a data.frame
#'
#' @export
#'
#' @aliases getTFBS
#' @rdname getTFBS

setMethod('getTFBS',  'TrenaValidator',

      function(obj, tbl.regions, fimo.threshold, conservation.threshold, pwmFile){
         tbl.match <- fimoBatch(tbl.regions, matchThreshold=fimo.threshold, genomeName="hg38",
                                pwmFile=pwmFile)
         tbl.match <- as.data.frame(gscores(phastCons7way.UCSC.hg38,
                                            GRanges(tbl.match)), stringsAsFactors=FALSE)
         tbl.match <- subset(tbl.match, default >= conservation.threshold)
         colnames(tbl.match)[1] <- "chrom"
         colnames(tbl.match)[grep("default", colnames(tbl.match))] <- "phast7"
         obj@state$regulatoryRegions <- tbl.match
         tbl.match
         })

#----------------------------------------------------------------------------------------------------
#' set the matrix to use in model building
#'
#' @param obj An instance of the TrenaValidator class
#' @param matrix a properly normalized numeric matrix with geneSymbol rownames
#'
#' @return noting
#'
#' @export
#'
#' @aliases setMatrix
#' @rdname setMatrix

setMethod('setMatrix',  'TrenaValidator',

      function(obj, matrix){
         obj@state$matrix <- matrix
      })

#----------------------------------------------------------------------------------------------------
#' return a trena model
#'
#' @param obj An instance of the TrenaValidator class
#'
#' @return a data.frame
#'
#' @export
#'
#' @aliases buildModel
#' @rdname buildModel

setMethod('buildModel',  'TrenaValidator',

      function(obj){
         tbl.reg <- obj@state$regulatoryRegions
         stopifnot(nrow(tbl.reg) > 0)
         solvers <- c("lasso", "ridge", "lassopv", "pearson", "spearman", "xgboost")
         x <- EnsembleSolver(obj@state$matrix,
                             obj@targetGene,
                             unique(tbl.reg$tf),
                             solverNames=solvers)

         tbl <- run(x)
         tbl.xtab <- as.data.frame(table(tbl.reg$tf))
         colnames(tbl.xtab) <- c("tf", "bindingSites")
         tbl <- merge(tbl, tbl.xtab, by.x="gene", by.y="tf")
         obj@state$model <- tbl
         tbl
         })

#----------------------------------------------------------------------------------------------------
