configfile: "config.yaml"

species = 'dm6'
bwa_index = config['knowledge']['indices'][species]['bwa']

data_root = '{person}/{project}/data/'

## read specific settings
reads_root = data_root + 'reads/'
read_suffix = '{sample}.fastq.gz'
raw_read_path = reads_root + 'raw/' 
raw_read_pattern= raw_read_path + read_suffix

mapping_root = data_root + 'mapping/'

# esat settings
esat_command = config['scripts']['esat_command']
esat_suffixes = ['window.txt', 'gene.txt']
esat_out_dir = data_root + 'esat/'
esat_out_prefix = esat_out_dir + '{sample}'
esat_out = [esat_out_prefix + '.' + suffix for suffix in esat_suffixes]

fastqc_suffixes = ['.zip', '.html']
fastqc_root = data_root + 'fastqc/'
fastqc_out = [fastqc_root + '{sample}_fastqc' + suffix for suffix in fastqc_suffixes]

m6a_motif = 'data/m6a_meme_motif.txt'
candidates_path = 'data/candidates/'
candidates_fa = candidates_path + 'candidates.fa'
control_seqs_fa = candidates_path + 'control_seqs.fa'

# find circ2
fc2 = config['scripts']['fc2']
find_circ2_root = mapping_root + 'circrna/' 
find_circ2_out_path = find_circ2_root + '{sample}/'
find_circ2_bwa_err = find_circ2_out_path + 'bwa.stderr.log'
find_circ2_bwa_pipe = find_circ2_out_path + 'bwa_mapped.pipe.bam'
find_circ2_err = find_circ2_out_path + 'fc.stderr.log'
find_circ2_run_dir = find_circ2_out_path + 'find_circ_run/'
find_circ2_out = [find_circ2_run_dir + suffix for suffix in ['circ_splice_sites.bed', 'lin_splice_sites.bed']]

# feature counts out
fc_out_path = data_root + 'counts/'
fc_out = fc_out_path + 'feature_counts.rds'

# include two snakefiles for star mapping and rfp read mapping
include: 'shared/code/snakefiles/mapWithStar.Snakefile'
include: 'shared/code/snakefiles/mapRFPReads.Snakefile'

### PREPARE FILES FOR DYNAMIC RULES ###
def get_out(person, project, to_expand):
	return expand(to_expand, person=person, project=project, sample=config['raw_data'][person][project]['samples'].keys())

def get_fc_in_pattern(wildcards):
	lib_type = config['raw_data'][wildcards.person][wildcards.project]['lib_type']
	if lib_type == 'rfp':
		pattern = rfp_mapped_out
	elif lib_type == 'star':
		pattern = star_out     
	else:	
		raise Exception('lib_type not valid for fc pattern')
	
	samples = config['raw_data'][wildcards.person][wildcards.project]['samples'].keys()

	return expand(pattern, person=wildcards.person, project=wildcards.project, sample=samples)

rule all:
	input:
		get_out('michela', 'm6a', fastqc_out),
		get_out('nagarjuna', 'rfp', fastqc_out),
		get_out('mor_osnat', 'synaptosomes', fastqc_out),
		get_out('michela', 'm6a', esat_out),
		get_out('michela', 'm6a', find_circ2_out),
		get_out('mor_osnat', 'synaptosomes', find_circ2_out),
		get_out('mor_osnat', 'synaptosomes', fc_out),
		get_out('nagarjuna', 'rfp', fc_out)

rule create_read_symlinks:
	output:
		raw_read_pattern
	params:
		p_link = lambda wildcards: config['raw_data'][wildcards.person][wildcards.project]['reads_dir'] + config['raw_data'][wildcards.person][wildcards.project]['samples'][wildcards.sample]
	shell:
		"ln -s {params.p_link} {output}"

rule create_gene_2_transcript_file:
    output:
        "shared/annotations/gene2transcript.txt"
    script:
        "shared/code/R/get_gene_2_transcript.R"

rule create_fastqc:
	input:
		raw_read_pattern,
	output:
		fastqc_out
	threads: 8
	params:
		fastqc_root = fastqc_root
	shell:
		"""
		mkdir -p {params.fastqc_root}

		fastqc -t {threads} -o {params.fastqc_root} {input}
		"""

rule index_bam_file:
	input:
		'{file}.bam'
	output:
		temp('{file}.bam.bai')
	shell:
		"samtools index {input}"

# rule get_read_length_dist:
# 	input:
# 		'{file}.bam'
# 	output:
# 		'{file}' + read_dist_suffix
# 	shell:
# 		 """
# 		 mkdir -p {filtered_out_path}
# 		 samtools view {input} | awk \'{{print(length($10))}}\' | sort | uniq -c > {output}
# 		 """

rule run_esat:
	input:
		annotation= rules.create_gene_2_transcript_file.output,
		star= star_out
	output:
		esat_out
	params:
		out_prefix= esat_out_prefix 
	shell:
		"""
		mkdir -p data/esat

		{esat_command} 	-in {input.star} \
			-geneMapping {input.annotation} \
			-out {params.out_prefix} \
			-task score3p
		"""
 
rule create_normalised_bigwig:
	input:
		mapped_reads=star_out,
		index=indexed_bam_out
	output:
		'data/bigwig/normalised/{sample}.{strand}.bw'
	shell:
		"""
		bamCoverage --bam {input.mapped_reads} -o {output} \
			--normalizeUsing RPKM \
			--binSize 5 \
			--filterRNAstrand {wildcards.strand} 
		"""

rule map_with_bwa_mem:
	input:
		raw_read_pattern
	output:
		bwa_out= pipe(find_circ2_bwa_pipe),
		bwa_err= find_circ2_bwa_err
	threads: 8
	shell:
		"""
		bwa mem -t {threads} -k 14 -T 1 -L 3,3 -O 6,6 -E 3,3 {bwa_index} {input} 2> {output.bwa_err} | \
		samtools view -Sbuh > {output.bwa_out}
		"""

rule find_circ2:
	input:
		find_circ2_bwa_pipe
	output:
		find_circ2_out
	threads: 8
	params:
		find_circ2_err= find_circ2_err,
		find_circ2_run_dir= find_circ2_run_dir
	conda:
		'shared/envs/find_circ2.yaml'
	shell:
		"""
		python2 {fc2} --bam --genome {bwa_index} --name {wildcards.sample} --output {params.find_circ2_run_dir} {input} 2> {params.find_circ2_err} 
		"""

rule count_features:
	input:
		get_fc_in_pattern
	output:
		fc_out		
	params:
		annotation= config['knowledge']['annotations'][species]
	threads: 16
	script:
		'shared/code/R/featureCounts.R'
