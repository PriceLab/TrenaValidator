library(TrenaValidator)
library(RUnit)
#------------------------------------------------------------------------------------------------------------------------
if(!exists("tv")) {
   message(sprintf("--- creating instance of TrenaValidator"))
   tbl.benchmark <- get(load(system.file(package="TrenaValidator", "extdata", "tbl.A.RData")))
   # mtx <- get(load(system.file(package="TrenaValidator", "extdata", "mtx.gtex.lung.RData")))
   tv <- TrenaValidator(TF="TWIST1", "MMP2", tbl.benchmark);
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

   tv <- TrenaValidator(TF="AHR", targetGene, tbl.benchmark);
   tbl.regions <- getSimplePromoter(tv)
   tbl.tfbs <- getTFBS(tv, tbl.regions, fimo.threshold=fimo, conservation.threshold=conservation, meme.file)
   tf.overlap <- unique(intersect(tbl.tfbs$tf, curated.tfs))
   printf("--- %d TF curated & predicted: %s", length(tf.overlap), paste(tf.overlap, collapse=", "))

   if(length(tf.overlap) == 0){
      printf("---- no binding sites include curated.tfs, abandoning %s at %5.2e, %5.2f", targetGene, fimo, conservation)
      return(NA)
      }


   #tbl.gh <- findEnhancers(tv, eliteOnly=TRUE)
   #dim(tbl.gh)
   #file.tfbs <- "tbl.tfbs.rela.RData"
   #if(file.exists(file.tfbs)){
   #   load(file.tfbs)
   #} else {
   #   save(tbl.tfbs, file=file.tfbs)
   #   }

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
         printf("--- model has %d predicted tfs", nrow(model))
         print(subset(tbl.model, gene %in% curated.tfs))
         break;   # one success is enough for now
         }
      } # for matrix.name

} # buildModelForTargetGene
#------------------------------------------------------------------------------------------------------------------------
tbl.targets <- as.data.frame(sort(table(tbl.benchmark$target)))
for(target in tbl.targets$Var1[6:10])
   buildModelForTargetGene(target, 1e-3, 0.75, 5000, 5000)
