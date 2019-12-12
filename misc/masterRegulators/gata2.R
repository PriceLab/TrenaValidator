library(TrenaValidator)
library(TrenaProjectHG38.generic)
library(BiocParallel)
library(RUnit)
library(rnaSeqNormalizer)
library(org.Hs.eg.db)
library(GO.db)
library(igvR)
library(GenomicScores)
library(phastCons7way.UCSC.hg38); phast.7 <- phastCons7way.UCSC.hg38
library (RColorBrewer)
library(PSICQUIC)
#------------------------------------------------------------------------------------------------------------------------
targetGene <- "GATA2"
#------------------------------------------------------------------------------------------------------------------------
if(!exists("tv")) {
   benchmark.full <- "~/github/trena/misc/saez-benchmark-paper/GarciaAlonso_Supplemental_Tables/database.csv"
   tbl.bm <-read.table(benchmark.full, sep=",", as.is=TRUE, header=TRUE, nrow=-1)

   message(sprintf("--- creating instance of TrenaValidator"))
   tbl.benchmark <- get(load(system.file(package="TrenaValidator", "extdata", "tbl.A.RData")))
   tbl.benchmark$pubmed.count <- unlist(lapply(strsplit(tbl.benchmark$pubmedID_from_curated_resources, ","), length))
   #mtx <- get(load(system.file(package="TrenaValidator", "extdata", "mtx.gtex.lung.RData")))
   tv <- TrenaValidator(TF="TWIST1", targetGene, tbl.benchmark);
   #setMatrix(tv, mtx)
   tp.hg38 <- TrenaProjectHG38.generic()
   }
#------------------------------------------------------------------------------------------------------------------------
if(!exists("igv")) {
   igv <- igvR()
   setGenome(igv, "hg38")
   setBrowserWindowTitle(igv, targetGene)
   }
#------------------------------------------------------------------------------------------------------------------------
# motifs <- query(MotifDb, "hsapiens", c("jaspar2018", "hocomoco"))
# meme.file <- "human.hocomoco.meme"
motifs <- query(MotifDb, "hsapiens", c("jaspar2018", "hocomoco"))
meme.file <- "human.jaspar2018.meme"
export(motifs, con=meme.file, format="meme")
#------------------------------------------------------------------------------------------------------------------------
if(!exists("genes.erythroid")){
   go.id <- "GO:0030218"  #  erythrocyte differentiation
   suppressWarnings(tbl.tfe <- select(org.Hs.eg.db, keys=go.id, keytype="GOALL", columns="SYMBOL"))
   genes.erythroid <- sort(unique(tbl.tfe$SYMBOL))
   length(genes.erythroid)
   }
if(!exists("genes.regulators")){
   genes.regulators <- c("GATA1", "GATA2", "ZFPM1", "KLF1", "FLI1", "TAL1", "CEBPA", "SPI1", "JUN",
                         "EGR1", "EGR2", "NAB1", "NAB2", "GFI1", "JUN", "HEY1")
   }

genes.other <- c("NFE2")
genes.all <- sort(unique(c(genes.regulators, genes.erythroid, genes.other)))
#------------------------------------------------------------------------------------------------------------------------
getCorcesMatrix <- function()
{
   file <- "~/github/TrenaProjectErythropoiesis/prep/import/buenrostro/GSE74246_RNAseq_All_Counts.txt"
   tbl <- read.table(file, sep="\t", as.is=TRUE, header=TRUE, nrow=-1)
   dim(tbl)
   rownames(tbl) <- tbl$X_TranscriptID
   mtx.counts <- as.matrix(tbl[, -1])
   fivenum(mtx.counts)
   normalizer <- rnaSeqNormalizer.gtex(mtx.counts, algorithm="log+scale", duplicate.selection.statistic="median")
   suppressWarnings(
      mtx.corces <- getNormalizedMatrix(normalizer)
      )
   x <- colnames(mtx.corces)
   mtx.corces[is.nan(mtx.corces)] <- 0
   deleters <- as.numeric(which(rowSums(abs(mtx.corces)) < 0.1))
   length(deleters)
   if(length(deleters) > 0)
      mtx.corces <- mtx.corces[-deleters,]
   fivenum(mtx.corces)  # -6.1187049 -0.5273388 -0.1581262  0.5399006  8.8888889
   dim(mtx.corces)      # 22438    81

   invisible(mtx.corces)

} # getCorcesMatrix
#------------------------------------------------------------------------------------------------------------------------
getBrandMatrix <- function()
{
   file <- "~/github/TrenaProjectErythropoiesis/prep/import/rnaFromMarjorie/GSE118537_DESeq_Read_Counts.tsv"
   tbl <- read.table(file, sep="\t", as.is=TRUE, header=TRUE, nrow=-1)
   dim(tbl)
   rownames(tbl) <- tbl$GeneName
   mtx.counts <- as.matrix(tbl[, -1])
   fivenum(mtx.counts)
   normalizer <- rnaSeqNormalizer.gtex(mtx.counts, algorithm="log+scale", duplicate.selection.statistic="median")
   suppressWarnings(
      mtx.brand <- getNormalizedMatrix(normalizer)
      )
   x <- colnames(mtx.brand)
   mtx.brand[is.nan(mtx.brand)] <- 0
   deleters <- as.numeric(which(rowSums(abs(mtx.brand)) == 0))
   length(deleters)
   if(length(deleters) > 0)
      mtx.brand <- mtx.brand[-deleters,]
   fivenum(mtx.brand)  # --4.82003077 -0.57054275 -0.06259517  0.62266976  5.10252039
   dim(mtx.brand)      # 23452    28

   invisible(mtx.brand)

} # getBrandMatrix
#------------------------------------------------------------------------------------------------------------------------
# run <- function(targetGene, phast7, fimo, upstream, downstream, display=FALSE)
# {
#
#    mtx <- getBrandMatrix()
#    tv <- TrenaValidator(TF=NA_character_, targetGene, tbl.benchmark);
#    setMatrix(tv, mtx)
#
#    #tbl.geneInfo <- getTranscriptsTable(tp.hg38, targetGene)
#    tbl.regions <- getSimplePromoter(tv, upstream=upstream, downstream=downstream)
#    browser()
#    showGenomicRegion(igv, with(tbl.regions, sprintf("%s:%d-%d", chrom, start, end)))
#    tbl.tfbs.fimo <- getTFBS.fimo(tv, tbl.regions, fimo.threshold=fimo, conservation.threshold=phast7, meme.file)
#    dim(tbl.tfbs.fimo)
#    tfs <- sort(unique(tbl.tfbs.fimo$tf))
#    tfs
#    length(tfs)
#
#    suppressWarnings(
#       tbl.model <- buildModel(tv, 0.1)
#       )
#    print(dim(tbl.model))
#    print(tbl.model)
#
#    tbl.tfbs <- subset(tbl.tfbs.fimo, tf %in% tbl.model$gene)
#
#    result <- list(model=tbl.model, tfbs=tbl.tfbs)
#    if(display)
#       displayModel(result)
#
#    result
#
# } # run
#------------------------------------------------------------------------------------------------------------------------
displayModel <- function(x)
{
   tbl.model <- x$model
   tbl.tfbs  <- x$tfbs
   shoulder <- 3000

   if("motif_id" %in% colnames(tbl.tfbs)){
      motifs <- unique(tbl.tfbs$motif_id)
      }
   if("motifName" %in% colnames(tbl.tfbs)){
      motifs <- unique(tbl.tfbs$motifName)
      }

   for(TF in tbl.model$gene){
     tbl.regions <- subset(tbl.tfbs, tf == TF)[, c("chrom", "start", "end", "tf")]
      track <- DataFrameAnnotationTrack(TF, tbl.regions, color="random", trackHeight=35)
      displayTrack(igv, track)
      } # for motif

} # displayModel
#------------------------------------------------------------------------------------------------------------------------
geneHancerTrack <- function(targetGene)
{
   setTargetGene(tp.hg38, targetGene)
   tbl.gh <- getEnhancers(tp.hg38)
   tbl.gh$width <- with(tbl.gh, 1 + end - start)
   #tbl.gh <- subset(tbl.gh, width < 5000)
   track <- DataFrameQuantitativeTrack("gh", tbl.gh[, c(1,2,3,11)], color="blue", autoscale=FALSE, min=0, max=50)
   displayTrack(igv, track)

}
#------------------------------------------------------------------------------------------------------------------------
conservationTrack <- function()
{
   loc <- getGenomicRegion(igv)
   starts <- with(loc, seq(start, end, by=5))
   ends <- starts + 5
   count <- length(starts)
   tbl.blocks <- data.frame(chrom=rep(loc$chrom, count), start=starts, end=ends, stringsAsFactors=FALSE)
   tbl.cons7 <- as.data.frame(gscores(phast.7, GRanges(tbl.blocks)), stringsAsFactors=FALSE)
   tbl.cons7$chrom <- as.character(tbl.cons7$seqnames)
   tbl.cons7 <- tbl.cons7[, c("chrom", "start", "end", "default")]
   track <- DataFrameQuantitativeTrack("phast7", tbl.cons7, autoscale=TRUE, color="red")
   displayTrack(igv, track)

} # conservationTrack
#------------------------------------------------------------------------------------------------------------------------
chipTrack <- function(tf, tissueRestriction=NA_character_, peaksAlso=FALSE)
{
   chrom.loc <- getGenomicRegion(igv)
   tbl.chip <- with(chrom.loc, getChipSeq(tp.hg38, chrom, start, end, tf))
   if(!is.na(tissueRestriction))
      tbl.chip <- tbl.chip[grep(tissueRestriction, tbl.chip$name),]

   printf("tf %s,  %d peaks", tf,  nrow(tbl.chip))

   if(nrow(tbl.chip) == 0){
      #printf("no ChIP for TF %s in %d bases", tf, chrom.loc$end - chrom.loc$start)
      return(data.frame())
      }
   tbl.track <- tbl.chip[, c("chrom", "start", "endpos", "name")]
   trackName <- sprintf("Ch-%s", tf)
   track <- DataFrameAnnotationTrack(trackName, tbl.track, color="random")
   displayTrack(igv, track)

   if(peaksAlso){
      tbl.track <- tbl.chip[, c("chrom", "peakStart", "peakEnd", "name")]
      trackName <- sprintf("peaks-%s", tf)
      track <- DataFrameAnnotationTrack(trackName, tbl.track, color="random")
      displayTrack(igv, track)
      } # peaksAlso

} # chipTrack
#------------------------------------------------------------------------------------------------------------------------
getATACseq <- function()
{
   roi <- getGenomicRegion(igv)
   chromosome <- roi$chrom
   start.loc <- roi$start
   end.loc <- roi$end

   directory <- "~/github/TrenaProjectErythropoiesis/prep/import/atacPeaks"
   files <- grep("narrowPeak$", list.files(directory), value=TRUE)
   result <- list()

   for(file in files){
      full.path <- file.path(directory, file)
      track.name <- sub("_hg38_macs2_.*$", "", sub("ATAC_Cord_", "", file))
      tbl.atac <- read.table(full.path, sep="\t", as.is=TRUE)
      colnames(tbl.atac) <- c("chrom", "start", "end", "name", "c5", "strand", "c7", "c8", "c9", "c10")
      tbl.atac.region <- subset(tbl.atac, chrom==chromosome & start >= start.loc & end <= end.loc)
      if(nrow(tbl.atac.region) > 0){
         tbl.atac.region$sample <- track.name
         result[[track.name]] <- tbl.atac.region
         }
      } # files

   tbl.out <- do.call(rbind, result)
   rownames(tbl.out) <- NULL

   tbl.out

} # getATACseq
#------------------------------------------------------------------------------------------------------------------------
displayATACseq <- function(union.only=FALSE)
{
   totalColorCount <- 12
   colors <- brewer.pal(8, "Dark2")
   currentColorNumber <- 0

   tbl.all <- getATACseq()
   samples <- unique(tbl.all$sample)
   current.day.string <- ""
   color <- colors[1]

   if(!union.only){
      for(current.sample in samples){
         this.day.string <- strsplit(current.sample, "_")[[1]][1]
         if(this.day.string != current.day.string){
            currentColorNumber <- (currentColorNumber %% totalColorCount) + 1
            color <- colors[currentColorNumber]
            current.day.string <- this.day.string
         }
         tbl.atac.sub <- subset(tbl.all, sample == current.sample)
         track.name <- current.sample
         track <- DataFrameQuantitativeTrack(track.name, tbl.atac.sub[, c("chrom", "start", "end", "c10")],
                                             color, autoscale=FALSE, min=0, max=430, trackHeight=30)
         displayTrack(igv, track)
         } # for samples
      } # if !union.only (that is: display each track)

   tbl.regions.condensed <- as.data.frame(union(GRanges(tbl.all[, c("chrom", "start", "end")]),
                                                GRanges(tbl.all[, c("chrom", "start", "end")])))[, c("seqnames", "start", "end")]
   colnames(tbl.regions.condensed) <- c("chrom", "start", "end")
   tbl.regions.condensed$chrom <- as.character(tbl.regions.condensed$chrom)
   lapply(tbl.regions.condensed, class)

   #state$tbl.regions.condensed <- tbl.regions.condensed
   track <- DataFrameAnnotationTrack("atac combined", tbl.regions.condensed, color="black")
   displayTrack(igv, track)

} # displayATACseq
#------------------------------------------------------------------------------------------------------------------------
run <- function(targetGene, phast7, fimo, bioc, upstream, downstream, display=FALSE)

{
   mtx <- getBrandMatrix()
   tv <- TrenaValidator(TF=NA_character_, targetGene, tbl.benchmark);
   setMatrix(tv, mtx)

   #tbl.geneInfo <- getTranscriptsTable(tp.hg38, targetGene)
   tbl.regions <- getSimplePromoter(tv, upstream=upstream, downstream=downstream)
   showGenomicRegion(igv, with(tbl.regions, sprintf("%s:%d-%d", chrom, start, end)))

   if(!is.na(fimo)){
      tbl.tfbs <- getTFBS.fimo(tv, tbl.regions, fimo.threshold=fimo, conservation.threshold=phast7, meme.file)
      print(dim(tbl.tfbs))
      tfs <- sort(unique(tbl.tfbs$tf))
      tfs
      length(tfs)
      }

   if(!is.na(bioc)){
      tbl.tfbs <- getTFBS.bioc(tv, tbl.regions, match.threshold=bioc, conservation.threshold=phast7, as.list(motifs))
      print(dim(tbl.tfbs))
      }

   suppressWarnings(
      tbl.model <- buildModel(tv, 1.0)
      )

   tbl.tfbs <- subset(tbl.tfbs, tf %in% tbl.model$gene)
   print(dim(tbl.model))
   print(tbl.model)

   tbl.model <- subset(tbl.model, abs(pearsonCoeff) > 0.4)
   result <- list(model=tbl.model, tfbs=tbl.tfbs)

   if(display)
      displayModel(result)

   result

} # run
#------------------------------------------------------------------------------------------------------------------------
tal1 <- function(conservation=0.5, fimo=NA_real_, bioc=NA_integer_, upstream=2500, downstream=500 )
{
  targetGene <- "TAL1"
  printf("modeling %s.  conservation: %5.2f   fimo: %5.2e   bioc: %d, upsteam:  %d   downstream: %d",
         targetGene, conservation, fimo, bioc, upstream, downstream)
  removeTracksByName(igv, getTrackNames(igv)[-1])
  x <- run(targetGene, conservation, fimo, bioc, upstream, downstream, TRUE)
  conservationTrack()
  displayATACseq(TRUE)
  junk <- lapply(x$model$gene, function(tf) chipTrack(tf, tissueRestriction="K562"))
  #geneHancerTrack(targetGene)
  invisible(x)

} # tal1
#------------------------------------------------------------------------------------------------------------------------
gata2 <- function(conservation=0.5, fimo=NA_real_, bioc=NA_integer_, upstream=2500, downstream=500 )
{
  targetGene <- "GATA2"
  printf("modeling %s.  conservation: %5.2f   fimo: %5.2e   bioc: %d, upsteam:  %d   downstream: %d",
         targetGene, conservation, fimo, bioc, upstream, downstream)
  tbl.regions <- getSimplePromoter(tv, upstream=upstream, downstream=downstream)
  showGenomicRegion(igv, with(tbl.regions, sprintf("%s:%d-%d", chrom, start, end)))

  removeTracksByName(igv, getTrackNames(igv)[-1])
  conservationTrack()
  x <- run(targetGene, conservation, fimo, bioc, upstream, downstream, TRUE)
  displayATACseq(FALSE)
  junk <- lapply(x$model$gene, function(tf) chipTrack(tf, tissueRestriction="K562"))
  geneHancerTrack(targetGene)
  for(tf in x$model$gene){
     chipTrack(tf, tissueRestriction="K562")
     }

  invisible(x)

} # gata2
#------------------------------------------------------------------------------------------------------------------------
explore.znf263.tfbs <- function()
{
   chrom <- "chr1"
   loc.start <- 47229390
   loc.end   <- 47229430
   seq <- as.character(getSeq(BSgenome.Hsapiens.UCSC.hg38, chrom, loc.start, loc.end))
    # "TGGAGAGAAA GAGGCAGGGCAAGAGGGAGGG AGAGAGAGAA"
                 "GGAGGAGGA?G?GGAGGAGG?"
    # consensusString(query(MotifDb, c("ZNF263", "jaspar2018"))[[1]])
    #        [GA]     [GA]       [AG]       G         [GC]        A          GAT       GAC        AG         GA
    #           1         2          3          4          5           6          7          8          9         10         11         12         13         14         15           16         17         18         19         20         21
    #A 0.15530030 0.3016738 0.72976698 0.02973416 0.06636035 0.913619954 0.10298654 0.10187069 0.53252379 0.33685592 0.20321628 0.43308172 0.27404004 0.26951099 0.52136528 2.509353e-01 0.20557926 0.54512635 0.34302593 0.26695110 0.47679685
    #C 0.01923203 0.0000000 0.00000000 0.00000000 0.11407942 0.036166721 0.05539875 0.09366590 0.06353791 0.10127995 0.07016738 0.03872662 0.04069577 0.10620282 0.06235642 6.563833e-05 0.00000000 0.00000000 0.01437479 0.04266492 0.05953397
    #G 0.72714145 0.6983262 0.18391861 0.90705612 0.81680341 0.002822448 0.73856252 0.75162455 0.36186413 0.46222514 0.62166065 0.44174598 0.58437808 0.59350181 0.31375123 6.340007e-01 0.75116508 0.38674106 0.58569084 0.64450279 0.40308500
    #T 0.09832622 0.0000000 0.08631441 0.06320971 0.00275681 0.047390876 0.10305218 0.05283886 0.04207417 0.09963899 0.10495569 0.08644568 0.10088612 0.03078438 0.10252708 1.149984e-01 0.04325566 0.06813259 0.05690843 0.04588119 0.06058418

   source("~/github/fimoService/batchMode/fimoBatchTools.R")
   meme.file <- "myTFS.meme"
   motifs <- query(MotifDb, c("jaspar2018", "ZNF263"))
   export(motifs, con=meme.file, format="meme")
   tbl.regions <- data.frame(chrom=chrom, start=loc.start, end=loc.end, stringsAsFactors=FALSE)

   tbl.match <- fimoBatch(tbl.regions, matchThreshold=1e-4, genomeName="hg38", pwmFile=meme.file)

    #   chrom    start      end     tf strand   score  p.value      matched_sequence                            motif_id
    # 1  chr1 47229397 47229417 ZNF263      + 15.5222 2.31e-06 GAAAGAGGCAGGGCAAGAGGG Hsapiens-jaspar2018-ZNF263-MA0528.1
    # 2  chr1 47229400 47229420 ZNF263      + 13.4667 8.47e-06 AGAGGCAGGGCAAGAGGGAGG Hsapiens-jaspar2018-ZNF263-MA0528.1
    # 3  chr1 47229401 47229421 ZNF263      + 17.2444 6.98e-07 GAGGCAGGGCAAGAGGGAGGG Hsapiens-jaspar2018-ZNF263-MA0528.1

    pfm <- motifs[[1]]
    Biostrings::maxScore(pfm)   # 13.67758
    hits <- Biostrings::matchPWM(pfm, seq, with.score=TRUE, min.score="79")
    hits
    mcols(hits)$score  # [1] 10.86518
    mcols(hits)$score/Biostrings::maxScore(pfm)  # 0.794


} # explore.znf263.tfbs
#------------------------------------------------------------------------------------------------------------------------
assessChIPPeak <- function(motif)
{
   tbl.regions <- as.data.frame(getGenomicRegion(igv), stringsAsFactors=FALSE)
   meme.file <- "tmp.meme"
   export(motif, con=meme.file, format="meme")
   tbl.fimo <- fimoBatch(tbl.regions, matchThreshold=1e-3, genomeName="hg38", pwmFile=meme.file)

   seq <- as.character(with(tbl.regions, getSeq(BSgenome.Hsapiens.UCSC.hg38, chrom, start, end)))
   hits <-  Biostrings::matchPWM(motif[[1]], seq, with.score=TRUE, min.score="50")
   tbl.bioc <- as.data.frame(ranges(hits))
   tbl.bioc$start <- tbl.bioc$start + tbl.regions$start
   tbl.bioc$end <- tbl.bioc$end + tbl.regions$start
   tbl.bioc$score <- mcols(hits)$score
   tbl.bioc$match <-  mcols(hits)$score/Biostrings::maxScore(motif[[1]])

   list(fimo=tbl.fimo, bioc=tbl.bioc)

} # assessChIPPeak
#------------------------------------------------------------------------------------------------------------------------
