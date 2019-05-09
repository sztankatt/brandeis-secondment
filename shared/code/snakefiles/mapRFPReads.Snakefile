missmatches = 0
# min and max readlength for mapped reads filtering
min_read_length = 27
max_read_length = 33

bowtie1_index = config['knowledge']['indices'][species]['bt1']
illumina_adaptors = config['knowledge']['adaptors']['illumina'] 
dm6_gtf = config['knowledge']['annotations'][species] 
dm6_tophat_tx_index = '/data/rajewsky/annotation/dm6/tophat_tx_index/ensGene.proper'
dm6_rRNA_bowtie2_index = config['knowledge']['indices'][species]['bt2_rRNA_plus']

# riboplus logs
ribo_plus_depleted_path = reads_root  + 'ribo_plus_depleted/'
ribo_plus_depleted_log = ribo_plus_depleted_path + '{sample}.log'
ribo_plus_depleted_pipe = ribo_plus_depleted_path + '{sample}.pipe.fastq.gz'
ribo_plus_depleted_mapped = ribo_plus_depleted_path + '{sample}.mapped.sam'

# tophat output
tophat_out_path = mapping_root + 'tophat/'
tophat_out = tophat_out_path + "{sample}"
tophat_out_mapped = tophat_out + '/accepted_hits.bam'
tophat_out_indexed = tophat_out_mapped + '.bai'

# trim short reads
trimmed_reads_path = reads_root + 'trimmed/'
trimmed_reads = trimmed_reads_path + read_suffix

# pattern of clean fastq.gz, trimmed and without umis
clean_reads_path = reads_root + 'clean/'
clean_reads = clean_reads_path + read_suffix

# extract 4 umis from both side
umi_extraction_regex = "'(?P<umi_1>.{4}).*(?P<umi_2>.{4})'"

deduped_out_path = mapping_root + 'deduped/'
deduped_out = deduped_out_path + '{sample}.bam'
deduped_log = deduped_out_path + '{sample}.log'

filtered_out_path = mapping_root + 'filtered/'
filtered_out = filtered_out_path + '{sample}_filtered.bam'

rfp_mapped_out = filtered_out

ruleorder:  filter_reads_length > remove_pcr_duplicates

rule trim_with_flexbar:
	input:
		raw_read_pattern
	output:
		pipe(trimmed_reads)
	threads: 4
	shell:
		"""
		mkdir -p {trimmed_reads_path}
		flexbar -qf i1.8 -qt 25 -m 20 --min-read-length 28 --adapter-trim-end RIGHT --adapter-preset SmallRNA \
			-n {threads} -r {input} -t {trimmed_reads_path}{wildcards.sample} -z GZ
		"""

rule get_umis:
	input:
		rules.trim_with_flexbar.output
	output:
		clean_reads
	params:
		log_file= lambda wildcards: clean_reads_path + wildcards.sample + ".umi_tools.log" 
	shell:
		"umi_tools extract --stdin={input} --bc-pattern={umi_extraction_regex} --extract-method=regex --log={params.log_file} --stdout={output}"
		
# deplete rRNA and other rnas. used for RFP
rule deplete_rRNA_rfp:
	input:
		rules.get_umis.output
	output:
		pipe_out= temp(ribo_plus_depleted_pipe),
		logfile= ribo_plus_depleted_log
	threads: 8
	shell:
		"bowtie2 -p {threads} -x {dm6_rRNA_bowtie2_index} -U {input} \
			--local -D 20 -R 3 -N 0 -L 16 -i S,1,0.5 \
			--un-gz {output.pipe_out} > /dev/null 2> {output.logfile}"

rule map_with_tophat:
	input:
		rules.deplete_rRNA_rfp.output.pipe_out
	output:
		temp(tophat_out_mapped)
	params:
		out_dir= lambda wildcards: expand(tophat_out,person =wildcards.person, project=wildcards.project, sample = wildcards.sample)
	threads: 8
	shell:
		"""
		tophat 	--no-coverage-search --bowtie1 --max-intron-length 260000 --num-threads {threads} --max-multihits 1 \
			--GTF {dm6_gtf} --no-novel-juncs --transcriptome-index {dm6_tophat_tx_index} --read-mismatches 0 \
			--library-type fr-firststrand --output-dir {params.out_dir} {bowtie1_index} {input}
		"""

rule remove_pcr_duplicates:
	input:
		mapped= tophat_out_mapped,
		index= tophat_out_indexed
	output:
		out= deduped_out,
		logfile= deduped_log
	shell:
		"umi_tools dedup --stdin={input.mapped} --log={output.logfile} > {output.out}"

rule filter_reads_length:
	input:
		deduped_out
	output:
		filtered_out
	shell:
		'samtools view -h {input} | awk \'length($10)>={min_read_length} && length($10)<={max_read_length} || $1 ~ /^@/ \' | samtools view -hb > {output}'

		
