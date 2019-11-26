file <- "run_matrix-brand_phast7-0.50_fimo-1.000000e-04_shoulder-5000"
if(!exist("x"))
    load(file)

#------------------------------------------------------------------------------------------------------------------------
add.targetGeneColumn <- function(model, targetGene)
{
  model$targetGene <- targetGene
  model
  
} # add.targetGeneColumn
#------------------------------------------------------------------------------------------------------------------------
models.all <- lapply(x, function(el) if("model" %in% names(el)) return (el$model) else return(data.frame()))
tbl.models <- do.call(rbind, models.all)
rownames(tbl.models) <- NULL
dim(tbl.models)

