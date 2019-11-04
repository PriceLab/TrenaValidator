#----------------------------------------------------------------------------------------------------
#' @import methods
#' @importFrom AnnotationDbi select
#' @import org.Hs.eg.db
#' @import TrenaProject
#' @import TrenaProjectErythropoiesis
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
                               matrix="matrix",
                               tbl.benchmark="data.frame",
                               state="environment",
                               trenaProject="TrenaProject",
                               quiet="logical"
                               ))

#----------------------------------------------------------------------------------------------------
setGeneric('findEnhancers',  signature='obj', function(obj, eliteOnly) standardGeneric('findEnhancers'))
setGeneric('getEnhancerTable',   signature='obj', function(obj) standardGeneric('getEnhancerTable'))
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

TrenaValidator <- function(TF, targetGene, mtx, tbl.benchmark, quiet=TRUE)

{
   state <- new.env(parent=emptyenv())
   tp <- TrenaProjectErythropoiesis()
   setTargetGene(tp, targetGene)
   tbl.geneInfo <- getTranscriptsTable(tp)
   tbl.fallbackEnhancers <- with(tbl.geneInfo,
                                 data.frame(chrom=chrom, start=tss-5000, end=tss+5000,
                                            stringsAsFactors=FALSE))
   state$enhancers <- tbl.fallbackEnhancers

   .TrenaValidator(TF=TF, targetGene=targetGene, matrix=mtx,
                   tbl.benchmark=tbl.benchmark, state=state, trenaProject=tp, quiet=quiet)

} # TrenaValidator, the constructor
#----------------------------------------------------------------------------------------------------
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
         obj@state$enhancers <- tbl.enhancers
         })

#----------------------------------------------------------------------------------------------------
#' return the enhancers
#'
#' @param obj An instance of the TrenaValidator class
#'
#' @return a data.frame
#'
#' @export
#'
#' @aliases getEnhancerTable
#' @rdname getEnhancerTable

setMethod('getEnhancerTable',  'TrenaValidator',

      function(obj){
         obj@state$enhancers
         })

#----------------------------------------------------------------------------------------------------
