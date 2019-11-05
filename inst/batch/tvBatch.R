library(TrenaValidator)
library(BiocParallel)
if(!exists("tbl.benchmark"))
   tbl.benchmark <- get(load(system.file(package="TrenaValidator", "extdata", "tbl.A.RData")))
mtx <- get(load(system.file(package="TrenaValidator", "extdata", "mtx.gtex.lung.RData")))

OUT_DIR <- "lung-eliteGH"
LOG_DIR <- "lung-eliteGH.log"

if(!file.exists(OUT_DIR))
   dir.create(OUT_DIR)

if(!file.exists(LOG_DIR))
   dir.create(LOG_DIR)
#------------------------------------------------------------------------------------------------------------------------
# steps:
#   get enhancers for targetGene
#   get high-scoring TFs in highly conserved regions of the enhancers
#   build trena model with those tfs and a GTEX expression matrix
#   save enhancers, motif matches, model
validateOnePair <- function(TF, targetGene)
{
   printf("--- validateOnePair(%s, %s)", TF, targetGene)

   tv <- TrenaValidator(TF, targetGene, mtx, tbl.benchmark)
   findEnhancers(tv, eliteOnly=TRUE)
   tbl.enhancers <- getEnhancerTable(tv)
   save(tbl.enhancers, file=sprintf("%s/%s-%s.RData", OUT_DIR, TF, targetGene))

} # validateOnePair
#------------------------------------------------------------------------------------------------------------------------
WORKERS <- 2
parallel <- TRUE

rows <- sample(seq_len(nrow(tbl.benchmark)),3)

if(parallel){
   bp.params <- MulticoreParam(stop.on.error=FALSE, log=TRUE, logdir=LOG_DIR, threshold="INFO", workers=WORKERS)
   results <- bptry({bplapply(rows,
                              function(r) with(tbl.benchmark[r,], validateOnePair(TF, target)),
                              BPPARAM=bp.params)})
  } else {
    results <- lapply(rows, function(r) with(tbl.benchmark[r,], validateOnePair(TF, target)))
    }


