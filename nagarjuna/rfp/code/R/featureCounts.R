library(Rsubread)

fc <- featureCounts(files = unlist(snakemake@input),
					annot.ext = snakemake@params$annotation,
					isGTFAnnotationFile = TRUE,
					strandSpecific = 1,
					nthreads = snakemake@threads)

saveRDS(fc, snakemake@output[[1]])
