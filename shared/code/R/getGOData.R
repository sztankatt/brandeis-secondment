library(topGO)
library(tidyverse)
library(magrittr)

getGOData <- function(genes_of_interest, control_genes, gene_id2go, ontology = 'MF'){
	# unify control and interest genes 
	all_genes <- c(genes_of_interest, control_genes) %>% unique

	# create 0s for interesgint genes and 1 for everything else
	gene_list <- factor(as.integer(all_genes %in% genes_of_interest))
	# name the elements of the list by geneid-s
	names(gene_list) <- all_genes

	new("topGOdata", ontology = ontology, allGenes = gene_list, annot = annFUN.gene2GO, gene2GO = gene_id2go,
		nodeSize = 10)
}
