/*
 * Workflow
 */

workflow COVERAGE {
	take:
	reads
	assembly // assembly needs to be prefixed with sample name if not co-assembly

	main:
	ref = BWA_INDEX(assembly)

	// only coassembly mode since we merge the assemblies in the other case
	reads_and_ref = reads.combine(ref).map{[it[0], it[1], it[3]]}
	aln = BWA_MEM(reads_and_ref)

	coverage = SAMTOOLS_COVERAGE(aln.bam.join(aln.bai))

	COMBINE_COVERAGE(coverage.collect{it[1]})

	emit:
	bam = BWA_MEM.out.bam
	bai = BWA_MEM.out.bai
	depth = COMBINE_COVERAGE.out.depth
}

/*
 * Processes
 */

process BWA_INDEX {
	tag "pooled"
	container "quay.io/biocontainers/bwa:0.7.3a--h7132678_7"
	label "medium_computation"

	input:
	tuple val(sample), path(assembly)

	output:
	tuple val(sample), path('ref*')

	script:
	"""
	bwa index ${assembly} -p ref
	"""
}

process BWA_MEM {
	tag "$sample"
	publishDir "$params.outdir/BWA", mode: 'copy'
	container "quay.io/biocontainers/mulled-v2-fe8faa35dbf6dc65a0f7f5d4ea12e31a79f73e40:c56a3aabc8d64e52d5b9da1e8ecec2031668596d-0"
	label "medium_computation"

	input:
	tuple val(sample), path(fastq), path(index)

	output:
	tuple val(sample), path('coverage*.bam'), emit: bam
	tuple val(sample), path('coverage*.bai'), emit: bai

	script:
	"""
	bwa mem -a -M -t $task.cpus ref $fastq \
		| samtools sort -@ $task.cpus -o coverage_${sample}.bam
	samtools index -@ $task.cpus coverage_${sample}.bam coverage_${sample}.bai
	"""
}

process SAMTOOLS_COVERAGE {
	tag "$sample"
	publishDir "$params.outdir/samtools", mode: 'copy'
	container "quay.io/biocontainers/samtools:1.15--h1170115_1"
	label "medium_computation"

	input:
	tuple val(sample), path(bam), path(bai)

	output:
	tuple val(sample), path("*.csv")

	script:
	"""
	samtools coverage --no-header -q $params.min_mapq $bam |
		awk 'BEGIN{FS="\\t";OFS=","}{print "$sample",\$1,\$6,\$7}' \\
		>> abund_${sample}.csv
	"""
}

process COMBINE_COVERAGE {
	tag "pooled"
	publishDir "$params.outdir/abundance_table", mode: 'copy'
	label "py_script"
	label "low_computation"

	input:
	path abund

	output:
	path "coverage.csv", emit: coverage
	path "depth.csv", emit: depth

	script:
	"""
	#!/usr/bin/env python

	import pandas as pd

	abund = pd.concat([
		pd.read_csv(f, header=None) for f in ['${abund.join("','")}']
	])
	abund.columns = ["sample", "contig", "coverage", "depth"]

	abund = abund.pivot("sample", "contig")

	abund.depth.to_csv("depth.csv")
	abund.coverage.to_csv("coverage.csv")
	"""
}
