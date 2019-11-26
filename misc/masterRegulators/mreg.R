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
#------------------------------------------------------------------------------------------------------------------------
if(!exists("tv")) {
   benchmark.full <- "~/github/trena/misc/saez-benchmark-paper/GarciaAlonso_Supplemental_Tables/database.csv"
   tbl.bm <-read.table(benchmark.full, sep=",", as.is=TRUE, header=TRUE, nrow=-1)

   message(sprintf("--- creating instance of TrenaValidator"))
   tbl.benchmark <- get(load(system.file(package="TrenaValidator", "extdata", "tbl.A.RData")))
   tbl.benchmark$pubmed.count <- unlist(lapply(strsplit(tbl.benchmark$pubmedID_from_curated_resources, ","), length))
   #mtx <- get(load(system.file(package="TrenaValidator", "extdata", "mtx.gtex.lung.RData")))
   tv <- TrenaValidator(TF="TWIST1", "MMP2", tbl.benchmark);
   #setMatrix(tv, mtx)
   tp.hg38 <- TrenaProjectHG38.generic()
   }
#------------------------------------------------------------------------------------------------------------------------
if(!exists("igv")) {
   igv <- igvR()
   setGenome(igv, "hg38")
   }
#------------------------------------------------------------------------------------------------------------------------
# motifs <- query(MotifDb, "hsapiens", c("jaspar2018", "hocomoco"))
# meme.file <- "human.hocomoco.meme"
motifs <- query(MotifDb, "hsapiens", c("jaspar2018"))
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
run <- function(targetGene, phast7, fimo, upstream, downstream, display=FALSE)

{
   mtx <- getBrandMatrix()
   tv <- TrenaValidator(TF=NA_character_, targetGene, tbl.benchmark);
   setMatrix(tv, mtx)

   #tbl.geneInfo <- getTranscriptsTable(tp.hg38, targetGene)
   tbl.regions <- getSimplePromoter(tv, upstream=2500, downstream=500)
   tbl.tfbs.fimo <- getTFBS.fimo(tv, tbl.regions, fimo.threshold=fimo, conservation.threshold=phast7, meme.file)
   dim(tbl.tfbs.fimo)
   tfs <- sort(unique(tbl.tfbs.fimo$tf))
   tfs
   length(tfs)

   #tbl.tfbs.bioc <- getTFBS.bioc(tv, tbl.regions, match.threshold=bioc.match, conservation.threshold=phast7, as.list(motifs))
   #dim(tbl.tfbs.bioc)
   #length(unique(tbl.tfbs.bioc$tf))

   suppressWarnings(
      tbl.model <- buildModel(tv, 0.1)
      )
   print(dim(tbl.model))
   print(tbl.model)

   tbl.tfbs <- subset(tbl.tfbs.fimo, tf %in% tbl.model$gene)

   result <- list(model=tbl.model, tfbs=tbl.tfbs)
   if(display)
      displayModel(result)

   result

} # run
#------------------------------------------------------------------------------------------------------------------------
displayModel <- function(x)
{
   tbl.model <- x$model
   tbl.tfbs  <- x$tfbs
   shoulder <- 3000

   showGenomicRegion(igv, sprintf("%s:%d-%d",
                                  tbl.tfbs$chrom[1],
                                  min(tbl.tfbs$start) - shoulder,
                                  max(tbl.tfbs$end) + shoulder))
   motifs <- unique(tbl.tfbs$motif_id)
   for(motif in motifs){
      tbl.regions <- subset(tbl.tfbs, motif_id == motif)[, c("chrom", "start", "end", "tf")]
      tf <- unique(tbl.regions$tf)
      track <- DataFrameAnnotationTrack(tf, tbl.regions, color="random", trackHeight=35)
      displayTrack(igv, track)
      } # for motif


   xyz <- 99

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
displayATACseq <- function()
{
   totalColorCount <- 12
   colors <- brewer.pal(8, "Dark2")
   currentColorNumber <- 0

   tbl.all <- getATACseq()
   samples <- unique(tbl.all$sample)
   current.day.string <- ""
   color <- colors[1]

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
run <- function(targetGene, phast7, fimo, upstream, downstream, display=FALSE)

{
   mtx <- getBrandMatrix()
   tv <- TrenaValidator(TF=NA_character_, targetGene, tbl.benchmark);
   setMatrix(tv, mtx)

   #tbl.geneInfo <- getTranscriptsTable(tp.hg38, targetGene)
   tbl.regions <- getSimplePromoter(tv, upstream=2500, downstream=500)
   tbl.tfbs.fimo <- getTFBS.fimo(tv, tbl.regions, fimo.threshold=fimo, conservation.threshold=phast7, meme.file)
   dim(tbl.tfbs.fimo)
   tfs <- sort(unique(tbl.tfbs.fimo$tf))
   tfs
   length(tfs)

   #tbl.tfbs.bioc <- getTFBS.bioc(tv, tbl.regions, match.threshold=bioc.match, conservation.threshold=phast7, as.list(motifs))
   #dim(tbl.tfbs.bioc)
   #length(unique(tbl.tfbs.bioc$tf))

   suppressWarnings(
      tbl.model <- buildModel(tv, 0.1)
      )
   print(dim(tbl.model))
   print(tbl.model)

   tbl.tfbs <- subset(tbl.tfbs.fimo, tf %in% tbl.model$gene)

   result <- list(model=tbl.model, tfbs=tbl.tfbs)
   if(display)
      displayModel(result)

   result

} # run
#------------------------------------------------------------------------------------------------------------------------
tal1 <- function()
{
  targetGene <- "TAL1"
  x <- run(targetGene, 0.9, 1e-5, 2500, 500, TRUE)
  conservationTrack()
  displayATACseq()
  geneHancerTrack(targetGene)
  invisible(x)

} # tal1
#------------------------------------------------------------------------------------------------------------------------
