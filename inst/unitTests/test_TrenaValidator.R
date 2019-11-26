library(TrenaValidator)
library(TrenaProjectHG38.generic)
library(BiocParallel)
library(RUnit)
library(rnaSeqNormalizer)
library(org.Hs.eg.db)
library(GO.db)
#------------------------------------------------------------------------------------------------------------------------
if(!exists("tv")) {
   benchmark.full <- "~/github/trena/misc/saez-benchmark-paper/GarciaAlonso_Supplemental_Tables/database.csv"
   tbl.bm <-read.table(benchmark.full, sep=",", as.is=TRUE, header=TRUE, nrow=-1)

   message(sprintf("--- creating instance of TrenaValidator"))
   tbl.benchmark <- get(load(system.file(package="TrenaValidator", "extdata", "tbl.A.RData")))
   tbl.benchmark$pubmed.count <- unlist(lapply(strsplit(tbl.benchmark$pubmedID_from_curated_resources, ","), length))
   mtx <- get(load(system.file(package="TrenaValidator", "extdata", "mtx.gtex.lung.RData")))
   tv <- TrenaValidator(TF="TWIST1", "MMP2", tbl.benchmark);
   setMatrix(tv, mtx)
   tp.hg38 <- TrenaProjectHG38.generic()
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
runTests <- function()
{
   test_constructor()
   test_enhancers()
   test_buildBindingSitesTable.fimo()
   test_buildBindingSitesTable.bioc()
   test_buildBindingSitesTable.bioc.parallel()
   test_buildModel()

} # runTests
#------------------------------------------------------------------------------------------------------------------------
test_constructor <- function()
{
   message(sprintf("--- test_constructor"))

   checkTrue(all(c("TrenaValidator") %in% is(tv)))

} # test_constructor
#------------------------------------------------------------------------------------------------------------------------
test_enhancers <- function()
{
   message(sprintf("--- test_enhancers"))
   tbl.gh <- findEnhancers(tv, eliteOnly=TRUE)
   checkTrue(nrow(tbl.gh) > 5)

} # test_constructor
#------------------------------------------------------------------------------------------------------------------------
test_buildBindingSitesTable.fimo <- function()
{
   message(sprintf("--- test_buildBindingSitesTable.fimo"))

   tbl.gh <- findEnhancers(tv, eliteOnly=TRUE)
   widths <- with(tbl.gh, end-start)
   keepers <- which(widths < 500)

   tbl.tfbs.1 <- getTFBS.fimo(tv, tbl.gh[1,], fimo.threshold=1e-3, conservation.threshold=0.95, meme.file)
   checkTrue(nrow(tbl.tfbs.1) < 5)

   tbl.tfbs.2 <- getTFBS.fimo(tv, tbl.gh[1,], fimo.threshold=1e-3, conservation.threshold=0, meme.file)
   checkTrue(nrow(tbl.tfbs.2) > 1000)

   tbl.tfbs.3 <- getTFBS.fimo(tv, tbl.gh[1,], fimo.threshold=1e-2, conservation.threshold=0.95, meme.file)
   checkTrue(nrow(tbl.tfbs.3) > 100)

   tbl.tfbs.4 <- getTFBS.fimo(tv, tbl.gh[1,], fimo.threshold=1e-3, conservation.threshold=0.75, meme.file)
   checkTrue(nrow(tbl.tfbs.4) > 5)

} # test_buildBindingSitesTable.fimo
#------------------------------------------------------------------------------------------------------------------------
test_buildBindingSitesTable.bioc <- function()
{
   message(sprintf("--- test_buildBindingSitesTable.bioc"))

   tbl.regions <- getSimplePromoter(tv, upstream=2500, downstream=500)
   tbl.tfbs.1 <- getTFBS.bioc(tv, tbl.regions, match.threshold=98, conservation.threshold=0.95, as.list(motifs))
   checkTrue(nrow(tbl.tfbs.1) < 10)

   tbl.tfbs.2 <- getTFBS.bioc(tv, tbl.regions, match.threshold=90, conservation.threshold=0.95, as.list(motifs))
   checkTrue(nrow(tbl.tfbs.2) > 50)

   checkTrue(all(tbl.tfbs.1$tf %in% tbl.tfbs.2$tf))

} # test_buildBindingSitesTable.bioc
#------------------------------------------------------------------------------------------------------------------------
test_buildBindingSitesTable.bioc.parallel <- function()
{
   message(sprintf("--- test_buildBindingSitesTable.bioc.parallel"))

   goi <- c("POLD2", "PTPRC", "SAMD4A")

   output.dir <- "~/github/trenaValidator/inst/unitTests/outDir"
   if(!file.exists(output.dir)) dir.create(output.dir)

   log.dir <- "~/github/trenaValidator/inst/unitTests/logDir"
   if(!file.exists(log.dir))
      dir.create(log.dir)


   f <- function(gene, upstream, downstream, match, conservation){
      tv <- TrenaValidator(TF="TWIST1", "POLD2", tbl.benchmark);
      tbl.regions <- getSimplePromoter(tv, upstream=upstream, downstream=downstream)
      tbl.tfbs <- getTFBS.bioc(tv, tbl.regions, match.threshold=match, conservation.threshold=conservation,
                               as.list(motifs))
      result <- list(gene=gene, upstream=upstream, downstream=downstream, match=match, conservation=conservation,
                     tfbs=tbl.tfbs)
      save.filename <- file.path(output.dir,
                                 sprintf("%s-%d.up-%d.down.RData", gene, upstream, downstream))
      save(result, file=save.filename)
      return(result)
      }

   time.lapply <- system.time(x <- lapply(goi, function(gene) f(gene, 2500, 2500, 95, 0.95)))

   bp.params <- MulticoreParam(stop.on.error=FALSE, log=TRUE, logdir=log.dir, threshold="INFO", workers=2)
   time.bplapply <- system.time(x.bp <- bplapply(goi, function(gene) f(gene, 2500, 2500, 95, 0.95),
                                                 BPPARAM=bp.params))

   checkTrue(as.list(time.bplapply)$elapsed < as.list(time.lapply)$elapsed)

} # test_buildBindingSitesTable.bioc.parallel
#------------------------------------------------------------------------------------------------------------------------
test_buildModel <- function()
{
   setMatrix(tv, mtx)

   tbl.regions <- getSimplePromoter(tv, upstream=2500, downstream=500)
   tbl.tfbs <- getTFBS.bioc(tv, tbl.regions, match.threshold=90, conservation.threshold=0.95, as.list(motifs))
   checkTrue("TWIST1" %in% tbl.tfbs$tf)

   dim(tbl.tfbs)
   setRegulatoryRegionsTable(tv, tbl.tfbs)

   tbl.model <- tryCatch({
      suppressWarnings(buildModel(tv))
      },
   error = function(e){
      printf(" failed model build: %s", e)
      return(data.frame())
      })

   printf("nrow(tbl.model): %d", nrow(tbl.model))
   checkTrue(nrow(tbl.model) > 8)
   checkTrue("TWIST1" %in% tbl.model$gene)

} # test_buildModel
#------------------------------------------------------------------------------------------------------------------------
test_CYP1A1 <- function()
{
   message(sprintf("--- test_CYP1A1"))
   curated.tfs <- subset(tbl.benchmark, target=="CYP1A1" & score == "A")$TF
   mtx <- get(load(system.file(package="TrenaValidator", "extdata", "mtx.gtex.lung.RData")))
   tv <- TrenaValidator(TF="AHR", "CYP1A1", tbl.benchmark);
   setMatrix(tv, mtx)
   tbl.gh <- findEnhancers(tv, eliteOnly=TRUE)
   tbl.tfbs <- getTFBS.fimo(tv, tbl.gh, fimo.threshold=1e-4, conservation.threshold=0.75, meme.file)
   suppressWarnings(
      tbl.model <- buildModel(tv)
      )
   dim(tbl.model)
   subset(tbl.model, gene %in% curated.tfs)

} # test_CYP1A1
#------------------------------------------------------------------------------------------------------------------------
test_anyTarget <- function()
{
   message(sprintf("--- test_anyTarget"))
   targetGene <- tbl.benchmark[sample(seq_len(nrow(tbl.benchmark)), 1), "target"]
   curated.tfs <- subset(tbl.benchmark, target==targetGene & score == "A")$TF
   printf("%s has %d curated tfs", targetGene, length(curated.tfs))
   mtx <- get(load(system.file(package="TrenaValidator", "extdata", "mtx.gtex.uterus.RData")))
   tv <- TrenaValidator(TF="AHR", targetGene, tbl.benchmark);
   setMatrix(tv, mtx)
   tbl.gh <- findEnhancers(tv, eliteOnly=TRUE)
   dim(tbl.gh)
   tbl.tfbs <- getTFBS.fimo(tv, tbl.gh, fimo.threshold=1e-3, conservation.threshold=0.25, meme.file)
   dim(tbl.tfbs)
   suppressWarnings(
      tbl.model <- buildModel(tv)
      )
   dim(tbl.model)
   subset(tbl.model, gene %in% curated.tfs)

} # test_anyTarget
#------------------------------------------------------------------------------------------------------------------------
test_NFE2 <- function()
{
   message(sprintf("--- test_NFE2"))
   targetGene <- "NFE2"
   targetGene <- "GATA1"
   #targetGene <- "RUNX1"
   #targetGene <- "E2F4"
   #mtx <- getCorcesMatrix()
   mtx <- getBrandMatrix()
   tv <- TrenaValidator(TF=NA_character_, targetGene, tbl.benchmark);
   setMatrix(tv, mtx)

   phast7 <- 0.90
   fimo <- 1e-4
   #bioc.match <- 95
   #shoulder <- 10000

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
      tbl.model <- buildModel(tv)
      )
   print(dim(tbl.model))
   print(tbl.model)

   tbl.model

} # test_nfe2
#------------------------------------------------------------------------------------------------------------------------
parameterized.run <- function(targetGene, matrix.name, phast7=0.75, fimo=1e-5, shoulder=5000)
{
   stopifnot(matrix.name %in% c("brand", "corces"))

   mtx <- switch(matrix.name,
                "brand" = getBrandMatrix(),
                "corces" = getCorcesMarix()
                )

   if(!targetGene %in% rownames(mtx)){
      return(list(gene=targetGene, model=data.frame(), bs=data.frame()))
      }

   tv <- TrenaValidator(TF="AHR", targetGene, tbl.benchmark);
   setMatrix(tv, mtx)

   tbl.geneInfo <- getTranscriptsTable(tp.hg38, targetGene)
   chrom <- tbl.geneInfo$chrom
   start.loc <- tbl.geneInfo$tss - shoulder
   end.loc <- tbl.geneInfo$tss + shoulder
   tbl.regions <- data.frame(chrom=chrom, start=start.loc, end=end.loc, end=end.loc, stringsAsFactors=FALSE)
   tbl.tfbs.fimo <- getTFBS.fimo(tv, tbl.regions, fimo.threshold=fimo, conservation.threshold=phast7, meme.file)
   dim(tbl.tfbs.fimo)
   length(unique(tbl.tfbs.fimo$tf))

   suppressWarnings(
      tbl.model <- buildModel(tv)
      )
   printf("---- %s", targetGene)
   tbl.model$targetGene <- targetGene
   tbl.model$phast7 <- phast7
   tbl.model$fimo <- fimo
   tbl.model$shoulder <- shoulder
   print(tbl.model)

   list(gene=targetGene, model=tbl.model, bs=tbl.tfbs.fimo)

} # parameterized.run
#------------------------------------------------------------------------------------------------------------------------
run <- function(goi=genes.all, phast7=0.5, fimo=1e-4, shoulder=5000, matrix="brand")
{

   f <- function(targetGene)
      tryCatch({
         parameterized.run(targetGene, matrix.name=matrix, phast7=phast7, fimo=fimo, shoulder=shoulder)
      },
      error=function(e){
         printf("problem with targetGene %s", targetGene)
         })

   x <- lapply(goi, f)
   names(x) <- goi
   filename <- sprintf("run_matrix-%s_phast7-%4.2f_fimo-%3e_shoulder-%d",
                       matrix, phast7, fimo, shoulder)

   save(x, file=filename)

} # run
#------------------------------------------------------------------------------------------------------------------------
if(!interactive()){
    #runTests()
    run()
    }
