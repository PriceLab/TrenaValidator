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
setGeneric('getSimplePromoter',  signature='obj', function(obj, upstream=5000, downstream=5000)
              standardGeneric('getSimplePromoter'))
setGeneric('findEnhancers',  signature='obj', function(obj, eliteOnly) standardGeneric('findEnhancers'))
setGeneric('getTFBS.fimo',  signature='obj', function(obj, tbl.regions, fimo.threshold, conservation.threshold, pwmFile)
              standardGeneric('getTFBS.fimo'))
setGeneric('getTFBS.bioc',  signature='obj', function(obj, tbl.regions, match.threshold, conservation.threshold, pwms)
              standardGeneric('getTFBS.bioc'))
setGeneric('getTFBS.moods',  signature='obj', function(obj, tbl.regions, match.threshold, conservation.threshold, motifs)
              standardGeneric('getTFBS.moods'))
setGeneric('setMatrix',   signature='obj', function(obj, matrix) standardGeneric('setMatrix'))
setGeneric('setRegulatoryRegionsTable',  signature='obj', function(obj, tbl.reg)
              standardGeneric('setRegulatoryRegionsTable'))

setGeneric('buildModel',  signature='obj', function(obj, geneCutoff=0) standardGeneric('buildModel'))
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
#' a simple single regulatory region oritened around the tss
#'
#' @param obj An instance of the TrenaValidator class
#' @param upstream numeric, default 5kb
#' @param downstream numeric, default 5kb
#'
#' @return data.frame
#'
#' @export
#'
#' @aliases getSimplePromoter
#' @rdname getSimplePromoter

setMethod('getSimplePromoter',  'TrenaValidator',

   function(obj, upstream=5000, downstream=5000) {
     tbl.geneInfo <- getTranscriptsTable(obj@trenaProject)
     tss <- tbl.geneInfo$tss
     start.loc <- tss - upstream
     end.loc   <- tss + downstream
     if(tbl.geneInfo$strand == -1){
        start.loc <- tss - downstream
        end.loc   <- tss + upstream
        }
     tbl.promoter <- with(tbl.geneInfo,
                          data.frame(chrom=chrom, start=start.loc, end=end.loc,
                                     stringsAsFactors=FALSE))
     tbl.promoter
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
            return(getSimplePromoter(obj))
         return(tbl.enhancers)
         })

#----------------------------------------------------------------------------------------------------
#' return the TFBS
#'
#' @param obj An instance of the TrenaValidator class
#' @param tbl.regions a data.frame with chrom, start, end columns
#' @param fimo.threshold numeric, 1e-4 for example
#' @param conservation.threshold numeric, between zero and one
#' @param pwmFile character string, a full or relative path to a file in meme format
#'
#' @return a data.frame
#'
#' @export
#'
#' @aliases getTFBS.fimo
#' @rdname getTFBS.fimo

setMethod('getTFBS.fimo',  'TrenaValidator',

      function(obj, tbl.regions, fimo.threshold, conservation.threshold, pwmFile){
         if(!obj@quiet) printf("starting TrenaValidator::getTFBS.fimo")
         tbl.match <- fimoBatch(tbl.regions, matchThreshold=fimo.threshold, genomeName="hg38",
                                pwmFile=pwmFile)
         if(!obj@quiet) printf("after TrenaValidator::getTFBS.fimo, fimoBatch")
         tbl.match <- as.data.frame(gscores(phastCons7way.UCSC.hg38,
                                            GRanges(tbl.match)), stringsAsFactors=FALSE)
         tbl.match <- subset(tbl.match, default >= conservation.threshold)
         colnames(tbl.match)[1] <- "chrom"
         colnames(tbl.match)[grep("default", colnames(tbl.match))] <- "phast7"
         obj@state$regulatoryRegions <- tbl.match
         tbl.match
         })

#----------------------------------------------------------------------------------------------------
#' return the TFBS
#'
#' @param obj An instance of the TrenaValidator class
#' @param tbl.regions a data.frame with chrom, start, end columns
#' @param match.threshold numeric, 1e-4 for example
#' @param conservation.threshold numeric, between zero and one
#' @param motifs  a list of motifs, the result of a MotifDb query
#'
#' @return a data.frame
#'
#' @export
#'
#' @aliases getTFBS.moods
#' @rdname getTFBS.moods

setMethod('getTFBS.moods',  'TrenaValidator',

      function(obj, tbl.regions, match.threshold, conservation.threshold, motifs){
         if(!obj@quiet) printf("starting TrenaValidator::getTFBS.moods")
         gr.regions <- with(tbl.regions, GRanges(seqnames=chrom, IRanges(start=start, end=end)))
         tbl.match <- matchMotif(MotifDb, motifs, "hg38", gr.regions, match.threshold, fimoDataFrameStyle=TRUE)
         if(!obj@quiet) printf("after TrenaValidator::getTFBS.moods, hits: %d", nrow(tbl.match))
         tbl.match <- as.data.frame(gscores(phastCons7way.UCSC.hg38,
                                            GRanges(tbl.match)), stringsAsFactors=FALSE)
         tbl.match <- subset(tbl.match, default >= conservation.threshold)
         colnames(tbl.match)[1] <- "chrom"
         colnames(tbl.match)[grep("default", colnames(tbl.match))] <- "phast7"
         tfs <- mcols(MotifDb[tbl.match$motif])$geneSymbol
         tbl.match$tf <- tfs
         obj@state$regulatoryRegions <- tbl.match
         tbl.match
         })

#----------------------------------------------------------------------------------------------------
#' return the TFBS using the bioc pwm match, trena MotifMatcher class
#'
#' @param obj An instance of the TrenaValidator class
#' @param tbl.regions a data.frame with chrom, start, end columns
#' @param match.threshold numeric, between 0 and 100
#' @param conservation.threshold numeric, between zero and one
#' @param pwms a list of motifs from MotifDb
#'
#' @return a data.frame
#'
#' @aliases getTFBS.bioc
#' @rdname getTFBS.bioc
#'
#' @export

setMethod('getTFBS.bioc',  'TrenaValidator',

      function(obj, tbl.regions, match.threshold, conservation.threshold, pwms){

         motifMatcher <- MotifMatcher(genomeName="hg38", pfms=pwms, quiet=obj@quiet)
         tbl.match <- findMatchesByChromosomalRegion(motifMatcher, tbl.regions,
                                                     pwmMatchMinimumAsPercentage=match.threshold)
         colnames(tbl.match)[grep("motifStart", colnames(tbl.match))] <- "start"
         colnames(tbl.match)[grep("motifEnd", colnames(tbl.match))] <- "end"

         tbl.match.scored <- as.data.frame(gscores(phastCons7way.UCSC.hg38,
                                                   GRanges(tbl.match)), stringsAsFactors=FALSE)
         tbl.match.scored <- subset(tbl.match.scored, default >= conservation.threshold)
         colnames(tbl.match.scored)[grep("seqnames", colnames(tbl.match.scored))] <- "chrom"
         colnames(tbl.match.scored)[grep("default", colnames(tbl.match.scored))] <- "phast7"
         tbl.match.scored$tf <- mcols(MotifDb[tbl.match.scored$motifName])$geneSymbol
         obj@state$regulatoryRegions <- tbl.match.scored
         tbl.match.scored
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
#' set the regulatory regions table
#'
#' @param obj An instance of the TrenaValidator class
#' @param tbl.reg  A data.frame
#'
#' @return noting
#'
#' @export
#'
#' @aliases setRegulatoryRegionsTable
#' @rdname setRegulatoryRegionsTable

setMethod('setRegulatoryRegionsTable',  'TrenaValidator',

      function(obj, tbl.reg){
         obj@state$regulatoryRegions <- tbl.reg
      })

#----------------------------------------------------------------------------------------------------
#' return a trena model
#'
#' @param obj An instance of the TrenaValidator class
#' @param geneCutoff numeric, betweeen 0 and 1, controls inclusion of weaker TFs in model
#'
#' @return a data.frame
#'
#' @export
#'
#' @aliases buildModel
#' @rdname buildModel

setMethod('buildModel',  'TrenaValidator',

      function(obj, geneCutoff=0){
         tbl.reg <- obj@state$regulatoryRegions
         stopifnot(nrow(tbl.reg) > 0)
         solvers <- c("lasso", "lassopv", "ridge", "randomforest", "pearson", "spearman", "xgboost")
         tfs <- unique(intersect(tbl.reg$tf, rownames(obj@state$matrix)))
         x <- EnsembleSolver(obj@state$matrix,
                             obj@targetGene,
                             tfs,
                             geneCutoff=geneCutoff,
                             solverNames=solvers)
         tbl <- run(x)
         tbl.xtab <- as.data.frame(table(tbl.reg$tf))
         colnames(tbl.xtab) <- c("tf", "bindingSites")
         tbl <- merge(tbl, tbl.xtab, by.x="gene", by.y="tf")
         if("lassoPValue" %in% colnames(tbl))
            tbl <- roundNumericColumnsInDataframe(tbl, 3, "lassoPValue")
         else
            tbl <- roundNumericColumnsInDataframe(tbl, 3, NA)
         obj@state$model <- tbl
         new.order <- order(abs(tbl$pearsonCoeff), decreasing=TRUE)
         tbl <- tbl[new.order,]
         rownames(tbl) <- NULL
         tbl
         })

#----------------------------------------------------------------------------------------------------
roundNumericColumnsInDataframe <- function(tbl, digits, pvalColumnNames=NA)
{
  tbl.pvals <- data.frame()
  tbl.main <- tbl

  if(!(all(is.na(pvalColumnNames)))){
     pval.cols <- grep(pvalColumnNames, colnames(tbl))
     stopifnot(length(pval.cols) == length(pvalColumnNames))
     tbl.pvals <- tbl[, pval.cols, drop=FALSE]
     tbl.main <- tbl[, -pval.cols, drop=FALSE]
     }

  numeric_columns <- sapply(tbl.main, mode) == 'numeric'
  tbl.main[numeric_columns] <-  round(tbl.main[numeric_columns], digits)

  tbl.out <- tbl.main
  if(ncol(tbl.pvals) > 0){
     tbl.pvals <- apply(tbl.pvals, 2, function(col) as.numeric(formatC(col, format = "e", digits = 2)))
     tbl.out <- cbind(tbl.out, tbl.pvals)[, colnames(tbl)]
     }

  tbl.out

} # roundNumericColumnsInDataframe
#------------------------------------------------------------------------------------------------------------------------
