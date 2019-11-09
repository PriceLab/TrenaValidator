library(TrenaValidator)
library(RUnit)
library(BiocParallel)
library(TrenaProjectHG38.generic)

LOGDIR <- "log.08nov19"
OUTDIR <- "out.09nov19"
if(!file.exists(OUTDIR))
   dir.create(OUTDIR)

if(!file.exists(LOGDIR))
   dir.create(LOGDIR)

#------------------------------------------------------------------------------------------------------------------------
tp <- TrenaProjectHG38.generic()
#------------------------------------------------------------------------------------------------------------------------
if(!exists("tv")) {
   message(sprintf("--- creating instance of TrenaValidator"))
   tbl.benchmark <- get(load(system.file(package="TrenaValidator", "extdata", "tbl.A.RData")))
   # mtx <- get(load(system.file(package="TrenaValidator", "extdata", "mtx.gtex.lung.RData")))
   tv <- TrenaValidator(TF="TWIST1", "MMP2", tbl.benchmark, quiet=FALSE);
   }
#------------------------------------------------------------------------------------------------------------------------
motifs <- query(MotifDb, "hsapiens", c("jaspar2018", "hocomoco"))
meme.file <- "human.hocomoco.meme"
export(motifs, con=meme.file, format="meme")
#------------------------------------------------------------------------------------------------------------------------
gtex.directory <- "~/github/GTEx/misc/prep"
tissue.matrix.files <- grep("RData$", list.files(gtex.directory), value=TRUE)
tissue.names <- sub(".RData", "", tissue.matrix.files)
matrix.list <- lapply(tissue.matrix.files, function(filename) {mtx <- get(load(file.path(gtex.directory, filename))); mtx})
length(tissue.matrix.files)
names(matrix.list) <- tissue.names
lapply(matrix.list, dim)
#------------------------------------------------------------------------------------------------------------------------
buildModelForTargetGene <- function(targetGene, fimo, conservation, upstream=5000, downstream=5000)
{
   curated.tfs <- subset(tbl.benchmark, target==targetGene & score == "A")$TF
   printf("%s has %d curated tf/s: %s", targetGene, length(curated.tfs), paste(curated.tfs, collapse=", "))


   recognizedGene <- setTargetGene(tp, targetGene)
   if(!recognizedGene){
       printf("--- no transcript info for %s, skipping", targetGene)
       return(NA)
       }

   tv <- TrenaValidator(TF="AHR", targetGene, tbl.benchmark, quiet=FALSE);
   tbl.regions <- getSimplePromoter(tv)
   printf("getting tfbs in %d regions", nrow(tbl.regions))
   tbl.tfbs <- getTFBS(tv, tbl.regions, fimo.threshold=fimo, conservation.threshold=conservation, meme.file)
   tf.overlap <- unique(intersect(tbl.tfbs$tf, curated.tfs))
   printf("--- %d TF curated & predicted: %s", length(tf.overlap), paste(tf.overlap, collapse=", "))

   if(length(tf.overlap) == 0){
      printf("---- no binding sites include curated.tfs, abandoning %s at %5.2e, %5.2f", targetGene, fimo, conservation)
      return(NA)
      }

   for(matrix.name in names(matrix.list)){
      printf("modeling %s with %s", targetGene, matrix.name)
      mtx <- matrix.list[[matrix.name]]
      setMatrix(tv, mtx)
      setRegulatoryRegionsTable(tv, tbl.tfbs)
      tbl.model <- tryCatch({
         tbl.model <- buildModel(tv)
         },
      error = function(e){
         printf(" failed model build: %s", e)
         return(data.frame())
         })
      if(any(curated.tfs %in% tbl.model$gene)){
         printf("--- model for %s has %d predicted tfs", targetGene, nrow(tbl.model))
         print(subset(tbl.model, gene %in% curated.tfs))
         result <- list(targetGene=targetGene, fimo=fimo, conservation=conservation,
                     upstream=upstream, downstream=downstream,
                     regions=tbl.tfbs,
                     model=tbl.model)
         results.file <- file.path(OUTDIR, sprintf("%s.RData", targetGene))
         printf("--- about to save results to %s", results.file)
         save(result, file=results.file)
         break;   # one success is enough for now
         }
      } # for matrix.name

} # buildModelForTargetGene
#------------------------------------------------------------------------------------------------------------------------
tbl.targets <- as.data.frame(sort(table(tbl.benchmark$target)))
parallel <- TRUE

#------------------------------------------------------------------------------------------------------------------------
run <- function()
{
   for(target in tbl.targets$Var1[207:210])
       buildModelForTargetGene(target, 1e-3, 0.75, 5000, 5000)

} # run
#------------------------------------------------------------------------------------------------------------------------
runParallel <- function()
{
   WORKERS <- 2
   #set.seed(17)
   rows <- sort(sample(seq_len(nrow(tbl.targets)), 100))
   #rows <- seq_len(nrow(tbl.targets))[300:350]
   parallel <- TRUE

   func <- function(row){
      targetGene <- as.character(tbl.targets$Var1[row])
      printf("func with targetGene: %s", targetGene)
      buildModelForTargetGene(targetGene, 1e-3, 0.75, 2500, 500)
      }

   if(parallel){
      bp.params <- MulticoreParam(stop.on.error=FALSE, log=TRUE, logdir=LOGDIR, threshold="INFO", workers=WORKERS)
      bp.params <- MulticoreParam(workers=WORKERS)
      #bp.params <- SerialParam()
      printf("starting bplapply")
      results <- bptry({bplapply(rows, func, BPPARAM=bp.params)})
      } # if parallel

   if(!parallel){
      results <- lapply(rows, func)
      }

} # runParallel
#------------------------------------------------------------------------------------------------------------------------
if(!interactive())
    runParallel()
