# function to read the output of star log
readStarLog <- function(log_file){
	out = list()
	lines = readLines(log_file)

	out$input_reads = (lines[6] %>% strsplit('\t') %>% unlist)[2] %>% as.integer

	out$uniq_reads = (lines[9] %>% strsplit('\t') %>% unlist)[2] %>% as.integer

	out$ratio = (lines[10] %>% strsplit('\t') %>% unlist)[2] %>% sub('%', '', .) %>% as.numeric
	out$ratio = out$ratio / 100

	out$avg_length = (lines[11] %>% strsplit('\t') %>% unlist)[2] %>% as.numeric
	
	return(out)
}
