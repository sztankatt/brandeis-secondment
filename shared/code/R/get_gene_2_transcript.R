####
# File which creates the gene 2 transcript mapping, required by the --geneMapping flag of esat short read counter
# TO BE CALLED BY SNAKEMAKE
####
library(biomaRt)
library(plyr)
library(tidyverse)

###
# TODO: add other species options, not only dm6
### 

ensembl = useMart('ensembl', dataset = 'dmelanogaster_gene_ensembl') 

# download all the exons, and rename the columns
exons <- ensembl %>%
	getBM(
		  attributes = c('ensembl_transcript_id',
						 'ensembl_gene_id',
						 'external_gene_name',
						 'exon_chrom_start', 
						 'exon_chrom_end',
						 'strand',
						 'chromosome_name',
						 'transcript_start', 
						 'transcript_end'),
		  mart = .) %>%
	as_tibble %>%
	rename(name = ensembl_transcript_id,
		   name2 = ensembl_gene_id,
		   chrom = chromosome_name,
		   exonStarts = exon_chrom_start,
		   exonEnds = exon_chrom_end,
		   txStart = transcript_start,
		   txEnd = transcript_end)

# mutate strand to -/+
exons <- exons %>%
	mutate(strand = ifelse(strand == 1, "+", "-"))

gene2transcript <- exons %>%
	# filter out 'artificial transcripts' where gene id and tx id is the same.
	# this was causing issues in esat
	filter(name != name2) %>%
	ddply(.,
		  .(name, name2, chrom, strand, txStart, txEnd),
		  summarize,
		  # collapse exon starts and ends, also add trailing comma
		  exonStarts = paste0(paste(exonStarts, collapse = ','), ','),
		  exonEnds = paste0(paste(exonEnds, collapse = ','), ',')) %>%
	as_tibble

write_tsv(gene2transcript, path = snakemake@output[[1]], col_names = T)
