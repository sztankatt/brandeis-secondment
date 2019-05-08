calculateEnrichment <- function(data, conditions){
# calculate ratios of enrichment of conditions against a given ccondition , default = input
	enrichment_data <- data %>%
		dplyr::select(condition, replicate, gene_id, gene_name, tpm) %>%
		# spread the tpm-s across conditions
		spread(condition, tpm) %>%
		# group tpm's back into one column, this time only IgG and IP
		gather('pull_down_condition', 'pull_down_tpm', conditions) %>%
		group_by(pull_down_condition, replicate) %>%
		mutate(tpm_enrichment = pull_down_tpm / input)

	# take back counts from inputs
	enrichment_data <- data %>%
		filter(condition == 'input') %>%
		ungroup() %>%
		transmute(gene_id = gene_id,
				  replicate = replicate,
				  input_count = count) %>%
		right_join(enrichment_data, by = c('gene_id', 'replicate'))

	enrichment_data
}
