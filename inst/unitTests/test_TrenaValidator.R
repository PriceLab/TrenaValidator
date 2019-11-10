library(TrenaValidator)
#library(TrenaProjectHG38.generic())
library(RUnit)
#------------------------------------------------------------------------------------------------------------------------
if(!exists("tv")) {
   message(sprintf("--- creating instance of TrenaValidator"))
   tbl.benchmark <- get(load(system.file(package="TrenaValidator", "extdata", "tbl.A.RData")))
   mtx <- get(load(system.file(package="TrenaValidator", "extdata", "mtx.gtex.lung.RData")))
   tv <- TrenaValidator(TF="TWIST1", "MMP2", tbl.benchmark);
   setMatrix(tv, mtx)
   }
#------------------------------------------------------------------------------------------------------------------------
motifs <- query(MotifDb, "hsapiens", c("jaspar2018", "hocomoco"))
meme.file <- "human.hocomoco.meme"
export(motifs, con=meme.file, format="meme")
#----------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_constructor()
   test_enhancers()
   test_buildBindingSitesTable.fimo()
   test_buildBindingSitesTable.bioc()
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

   tbl.tfbs.1 <- getTFBS(tv, tbl.gh[1,], fimo.threshold=1e-3, conservation.threshold=0.95, meme.file)
   checkTrue(nrow(tbl.tfbs.1) < 5)

   tbl.tfbs.2 <- getTFBS(tv, tbl.gh[1,], fimo.threshold=1e-3, conservation.threshold=0, meme.file)
   checkTrue(nrow(tbl.tfbs.2) > 1000)

   tbl.tfbs.3 <- getTFBS(tv, tbl.gh[1,], fimo.threshold=1e-2, conservation.threshold=0.95, meme.file)
   checkTrue(nrow(tbl.tfbs.3) > 100)

   tbl.tfbs.4 <- getTFBS(tv, tbl.gh[1,], fimo.threshold=1e-3, conservation.threshold=0.75, meme.file)
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
   message(sprintf("--- test_buildBindingSitesTable.bioc"))

   goi <- c("POLD2", "PTPRC", "SAMD4A")

   f <- function(gene, upstream, downstream, match, conservation){
      tv <- TrenaValidator(TF="TWIST1", "POLD2", tbl.benchmark);
      tbl.regions <- getSimplePromoter(tv, upstream=upstream, downstream=downstream)
      tbl.tfbs <- getTFBS.bioc(tv, tbl.regions, match.threshold=match, conservation.threshold=conservation, as.list(motifs))
      result <- list(gene=gene, upstream=upstream, downstream=downstream, match=match, conservation=conservation, tfbs=tbl.tfbs)
      save.filename <- file.path("~/github/trenaValidator/inst/unitTests/outDir",
                                 sprintf("%s-%d.up-%d.down.RData", gene, upstream, downstream))
      save(result, file=save.filename)
      return(result)
      }

   time.lapply <- system.time(x <- lapply(goi, function(gene) f(gene, 2500, 2500, 95, 0.95)))
   log.dir <- "~/github/trenaValidator/inst/unitTests/logDir"
   bp.params <- MulticoreParam(stop.on.error=FALSE, log=TRUE, logdir=log.dir, threshold="INFO", workers=2)
   time.bplapply <- system.time(x.bp <- bplapply(goi, function(gene) f(gene, 2500, 2500, 95, 0.95), BPPARAM=bp.params))

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
   tbl.tfbs <- getTFBS(tv, tbl.gh, fimo.threshold=1e-4, conservation.threshold=0.75, meme.file)
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
   tbl.tfbs <- getTFBS(tv, tbl.gh, fimo.threshold=1e-3, conservation.threshold=0.25, meme.file)
   dim(tbl.tfbs)
   suppressWarnings(
      tbl.model <- buildModel(tv)
      )
   dim(tbl.model)
   subset(tbl.model, gene %in% curated.tfs)

} # test_anyTarget
#------------------------------------------------------------------------------------------------------------------------
if(!interactive())
   runTests()
