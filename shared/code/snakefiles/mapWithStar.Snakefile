bt2_rRNA_index = config['knowledge']['indices'][species]['bt2_rRNA']

# ribo depletion
ribo_depleted_root= reads_root + 'ribo_depleted/'
ribo_depleted_out= ribo_depleted_root + '{sample}.depleted.fastq.gz'
ribo_depleted_log= ribo_depleted_root + '{sample}.log'

# star settings
star_suffix = 'Aligned.sortedByCoord.out.bam'
star_out_dir = mapping_root + 'star/'
star_out_prefix = star_out_dir + '{sample}.'
star_out = star_out_prefix + star_suffix
star_err_file = star_out_prefix + 'stderr.txt'
indexed_bam_out = star_out + '.bai'

rule deplete_rRNA:
	input:
		raw_read_pattern
	output:
		reads= temp(ribo_depleted_out),
		logfile= ribo_depleted_log
	threads: 8
	shell:
		"""
		bowtie2 -p {threads} -x {bt2_rRNA_index} -U {input} \
			--very-fast-local \
			--un-gz {output.reads} > /dev/null 2> {output.logfile}
		"""

rule map_with_star:
	input:
		reads= rules.deplete_rRNA.output.reads,
		index= config['knowledge']['indices'][species]['star']
	output:
		star_out
	threads: 8
	params:
		mapping_dir= star_out_dir,
		out_prefix= star_out_prefix,
		err_file= star_err_file
	shell:
		"""
		mkdir -p {params.mapping_dir}

		STAR    --runThreadN {threads} \
				--outFileNamePrefix {params.out_prefix} \
				--runMode alignReads \
				--genomeDir {input.index} \
				--readFilesIn {input.reads} \
				--readFilesCommand zcat \
				--outSAMtype BAM SortedByCoordinate	2>{params.err_file}
		"""
