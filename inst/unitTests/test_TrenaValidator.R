library(TrenaValidator)
library(RUnit)
#------------------------------------------------------------------------------------------------------------------------
if(!exists("tProj")) {
   message(sprintf("--- creating instance of TrenaValidator"))
   tbl.benchmark <- get(load(system.file(package="TrenaValidator", "extdata", "tbl.A.RData")))
   mtx <- get(load(system.file(package="TrenaValidator", "extdata", "mtx.gtex.lung.RData")))
   tv <- TrenaValidator(TF="TWIST1", "MMP2", mtx, tbl.benchmark);
   }
#------------------------------------------------------------------------------------------------------------------------
runTests <- function()
{
   test_constructor()

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
   tbl.fallback.enhancers <- getEnhancerTable(tv)
   checkEquals(with(getEnhancerTable(tv), end-start), 10000)

   findEnhancers(tv, eliteOnly=TRUE)
   tbl.enhancers <- getEnhancerTable(tv)
   checkTrue(nrow(tbl.enhancers) > 9)
   checkTrue(sum(with(tbl.enhancers, end-start)) > 25000)

} # test_constructor
#------------------------------------------------------------------------------------------------------------------------
if(!interactive())
   runTests()
