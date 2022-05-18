include { READS_QC as QC_BEFORE } from './reads_qc.nf' \
	addParams(outdir: "${params.outdir}/qc_before")
include { READS_QC as QC_AFTER } from './reads_qc.nf' \
	addParams(outdir: "${params.outdir}/qc_after")

/*
 * Workflows
 */

workflow TRIMMING {
	take: reads

	main:

	QC_BEFORE(reads)

	trimmed = params.trimmer.toLowerCase() == "fastp" ?
		FASTP(reads) :
		TRIMMOMATIC(reads)

	QC_AFTER(trimmed.paired.mix(trimmed.unpaired))

	emit:
	paired = trimmed.paired
	unpaired = trimmed.unpaired
}

/*
 * Processes
 */

process TRIMMOMATIC {
	tag "$sample"
	publishDir "${params.outdir}/Trimmomatic", mode: 'copy'
	container "quay.io/biocontainers/trimmomatic:0.32--hdfd78af_4"
	label "medium_computation"

	input:
	tuple val(sample), path(fastqs)

	output:
	tuple val(sample), path("paired/*_R{1,2}.fastq.gz"), emit: paired
	tuple val(sample), path("unpaired/*_unpaired.fastq.gz"), optional: true, emit: unpaired

	script:
	adapter_arg = params.adapters == "" ? "" : "ILLUMINACLIP:$params.adapters:2:30:10:2:True"
	cmd = fastqs[0].getExtension() == "gz" ? "zcat" : "cat"
	"""
	#!/usr/bin/env bash
	read_len=`$cmd ${fastqs[0]} | head -2 | tail -1 | wc -c`
	crop_tail=\$((\$read_len-$params.trim_tail1))

	mkdir paired unpaired

	trimmomatic PE -threads $task.cpus $fastqs \\
	  paired/${sample}_R1.fastq.gz unpaired/${sample}.fastq.gz \\
	  paired/${sample}_R2.fastq.gz unpaired/${sample}.fastq.gz \\
	  $adapter_arg HEADCROP:$params.trim_front1 CROP:\$crop_tail \\
	  SLIDINGWINDOW:$params.cut_window_size:$params.cut_mean_quality
	"""
}

process FASTP {
	tag "$sample"
	publishDir "${params.outdir}/FASTP", mode: 'copy'
	container "quay.io/biocontainers/fastp:0.23.2--hb7a2d85_2"
	label "medium_computation"

	input:
	tuple val(sample), path(fastqs)

	output:
	tuple val(sample), path("paired/*_R{1,2}.fastq.gz"), emit: paired
	tuple val(sample), path("unpaired/*_unpaired.fastq.gz"), optional: true, emit: unpaired
	path "*.html", emit: html

	script:
	adapter_arg = params.adapters == "" ?
		"--detect_adapter_for_pe" :
		"--adapter_fasta $params.adapters"
	"""
	#!/usr/bin/env bash

	mkdir paired unpaired

	fastp --correction -w $task.cpus $adapter_arg --cut_right \\
	  --cut_window_size ${params.cut_window_size} --cut_mean_quality ${params.cut_mean_quality} \\
	  --trim_front1 ${params.trim_front1} --trim_front2 ${params.trim_front2} \\
	  --trim_tail1 ${params.trim_tail1} --trim_tail2 ${params.trim_tail2} \\
	  -i ${fastqs[0]} -I ${fastqs[1]} \\
	  -o paired/${sample}_R1.fastq.gz -O paired/${sample}_R2.fastq.gz \\
	  --unpaired1 unpaired/${sample}_unpaired.fastq.gz \\
	  --unpaired2 unpaired/${sample}_unpaired.fastq.gz \\
	  --html fastp_${sample}.html

	if [ ! -s unpaired/${sample}_unpaired.fastq.gz ]; then
		rm -f unpaired/${sample}_unpaired.fastq.gz
	fi
	"""
}
