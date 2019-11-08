library(TrenaValidator)
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
   test_buildBindingSitesTable()
   #test_buildModel()

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
test_buildBindingSitesTable <- function()
{
   message(sprintf("--- test_buildBindingSitesTable"))

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

} # test_buildBindingSitesTable
#------------------------------------------------------------------------------------------------------------------------
test_buildModel <- function()
{
   message(sprintf(printf("--- test_buildModel")))
   tbl.gh <- findEnhancers(tv, eliteOnly=TRUE)
   tbl.tfbs.1 <- getTFBS(tv, tbl.gh, fimo.threshold=1e-5, conservation.threshold=0.95, meme.file)

   dim(tbl.tfbs.1)
   length(unique(tbl.tfbs.1$tf))

   suppressWarnings(
      tbl.model.1 <- buildModel(tv)
      )
   dim(tbl.model.1)

   tbl.tfbs.2 <- getTFBS(tv, tbl.gh, fimo.threshold=1e-3, conservation.threshold=0, meme.file)
   dim(tbl.tfbs.2)
   length(unique(tbl.tfbs.2$tf))

   suppressWarnings(
      tbl.model.2 <- buildModel(tv)
      )
   dim(tbl.model.2)

   tbl.reg <- tv@state$regulatoryRegions
   dim(tbl.reg)
   tf.keepers <- unique(subset(tbl.reg, phast7 > 0.9 & p.value < 1e-5)$tf)
   length(tf.keepers)
   tbl.model.2.refiltered <- subset(tbl.model.2, gene %in% tf.keepers)

   tbl.reg.filtered <- dim(subset(tv@state$regulatoryRegions, phast7 > 0.9 && p.value < 1e-4))

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
