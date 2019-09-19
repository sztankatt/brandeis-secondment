library(readr)
library(tidyverse)

parseFimoOut <- function(fimo_table, feature_lengths, tx2gene, nbins = 30){
	motifs <- fimo_table %>% 
		mutate(start = as.integer(start))
	
	motif_coverage <- motifs %>%
		dplyr::select(tx_id, start, motif_id) %>%
		inner_join(tx2gene, by='tx_id') %>%
		inner_join(feature_lengths, by = 'tx_id') %>%
		mutate(fivep_start = 0,
			   cds_start = fivep,
			   threep_start = fivep + cds)

	fivep_coverage <- motif_coverage %>%
		filter(start < cds_start) %>%
		mutate(start = start / fivep - 1) %>%
		gather('feature_type', 'feature_length', threep:cds) %>%
		filter(feature_type == 'fivep')

	cds_coverage <- motif_coverage %>%
		filter(start >= cds_start, start < threep_start) %>%
		mutate(start = (start - cds_start) / cds) %>%
		gather('feature_type', 'feature_length', threep:cds) %>%
		filter(feature_type == 'cds')

	threep_coverage <- motif_coverage %>%
		filter(start > threep_start) %>%
		mutate(start = 1 + (start - threep_start) / threep) %>%
		gather('feature_type', 'feature_length', threep:cds) %>%
		filter(feature_type == 'threep')

	all_coverage <- bind_rows(fivep_coverage,
			  cds_coverage,
			  threep_coverage)

	all_coverage <- all_coverage %>%
		mutate(bin = floor((start + 1) * nbins / 3) + 0.5) %>%
		group_by(gene_name, bin) %>% 
		mutate(norm_count = n()) %>%
		group_by(gene_name) %>%
		# normalise the count by the number of motifs in a given gene
		mutate(norm_count = norm_count  / (sum(norm_count)))

	all_coverage
}
