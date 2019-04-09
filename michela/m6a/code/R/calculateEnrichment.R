calculateEnrichment <- function(data, conditions, against = 'input'){
# calculate ratios of enrichment of conditions against a given ccondition , default = input
	enrichment_data <- data %>%
		dplyr::select(condition, replicate, gene_id, gene_name, tpm) %>%
		# spread the tpm-s across conditions
		spread(condition, tpm) %>%
		# group tpm's back into one column, this time only IgG and IP
		gather('pull_down_condition', 'pull_down_tpm', conditions) %>%
		group_by(pull_down_condition, replicate) %>%
		mutate(tpm_enrichment = pull_down_tpm / input)

	enrichment_data
}
