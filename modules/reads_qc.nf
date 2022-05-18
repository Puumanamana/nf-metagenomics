/*
 * Workflows
 */

workflow READS_QC {
	take:
	reads

	main:
	FASTQC(reads)
	MULTIQC(FASTQC.out.zip.collect())
}

/*
 * Processes
 */

process FASTQC {
	tag "$sample"
	publishDir "${params.outdir}/FASTQC", pattern: "*.html", mode: 'copy'
	container 'quay.io/biocontainers/fastqc:0.11.9--hdfd78af_1'
	label "low_computation"

	input:
	tuple val(sample), path(fastq)

	output:
	path "*.zip", emit: zip
	path "*.html", emit: html

	script:
	"""
	fastqc -t ${task.cpus} *.fastq*
	"""
}

process MULTIQC {
	tag "all"
	publishDir "${params.outdir}/MULTIQC", pattern: "*.html", mode: 'copy'
	container 'quay.io/biocontainers/multiqc:1.12--pyhdfd78af_0'
	label "low_computation"

	input:
	path fastqs

	output:
	path "multiqc.html"

	script:
	"""
	multiqc . -n multiqc.html
	"""
}
