library(TrenaValidator)
library(RUnit)
#------------------------------------------------------------------------------------------------------------------------
if(!exists("tv")) {
   message(sprintf("--- creating instance of TrenaValidator"))
   tbl.benchmark <- get(load(system.file(package="TrenaValidator", "extdata", "tbl.A.RData")))
   mtx <- get(load(system.file(package="TrenaValidator", "extdata", "mtx.gtex.lung.RData")))
   tv <- TrenaValidator(TF="TWIST1", "MMP2", mtx, tbl.benchmark);
   }
#------------------------------------------------------------------------------------------------------------------------
motifs <- query(MotifDb, c("hsapie"), c("jaspar2018", "hocomoco"))
export(motifs, con="jaspar2018-hocomoc-human.meme", format="meme")
#----------------------------------------------------------------------------------------------------

runTests <- function()
{
   test_constructor()
   test_enhancers()
   test_buildBindingSitesTable()

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
   message(sprintf(printf("--- test_buildBindingSitesTable")))
   motifs <- query(MotifDb, "hsapiens", c("jaspar2018", "hocomoco"))
   meme.file <- "human.hsapiens.hocomoco.meme"
   export(motifs, con=meme.file, format="meme")

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
if(!interactive())
   runTests()
