include { SEQTK_SEQ as MIN_LENGTH_FILTER } from './seqtk.nf'

/*
 * Workflow
 */

workflow ASSEMBLY {
	take:
	reads

	main:
	if(params.coassembly) {
		// concatenate forward and reverse fastqs
		reads = reads
			.transpose()
			.collectFile(sort: true) { it ->
				def suffix = (it[1].getName() =~/_([^_]*.fastq(.gz)?)/)[0][1];
				["all_$suffix", it[1]]
			}.map{['pooled', it]}.groupTuple()
	} else {
		reads = reads.transpose().groupTuple()
	}

	assembly = params.assembler == 'megahit' ?
		MEGAHIT(reads) :
		SPADES(reads)

	QUAST(assembly.collect{it[1]})

	// Length filtering
	assembly_filt = MIN_LENGTH_FILTER(assembly)

	emit:
	raw=assembly
	filt=assembly_filt
}

/*
 * Processes
 */

process SPADES {
	tag "$sample"
	publishDir "${params.outdir}/${params.assembler}", mode: 'copy'
	container "quay.io/biocontainers/spades:3.15.4--h95f258a_0"
	label "high_computation"

	input:
	tuple val(sample), path(fastq)

	output:
	tuple val(sample), path("assembly_*.fasta")

	script:
	"""
	#!/usr/bin/env bash

	inputs="-1 \$(ls *_R1.fastq.gz) -2 \$(ls *_R2.fastq.gz)"

	if compgen -G "*_unpaired.fastq.gz" > /dev/null; then
		inputs+=" -s \$(ls *_unpaired.fastq.gz)"
	fi

	${params.assembler}.py \$inputs -t $task.cpus -m ${task.memory.getGiga()} -o ${params.assembler}_output

	mv ${params.assembler}_output/scaffolds.fasta "assembly_${params.assembler}_${sample}.fasta"
	"""
	}

process MEGAHIT {
	tag "$sample"
	publishDir "${params.outdir}/megahit", mode: 'copy'
	container "quay.io/biocontainers/megahit:1.2.9--h2e03b76_1"
	label "high_computation"

	input:
	tuple val(sample), path(fastq)

	output:
	tuple val(sample), path("assembly_*.fasta")

	script:
	"""
	#!/usr/bin/env bash

	inputs="-1 \$(ls *_R1.fastq.gz) -2 \$(ls *_R2.fastq.gz)"
	if compgen -G "*_unpaired.fastq.gz" > /dev/null; then
		inputs+=" -r \$(ls *_unpaired.fastq.gz)"
	fi

	megahit -o megahit_out \$inputs

	mv megahit_out/final.contigs.fa assembly_megahit_${sample}.fasta
	"""
}

process QUAST {
	tag "all"
	publishDir "${params.outdir}/quast", mode: 'copy'
	errorStrategy "ignore"
	container "quay.io/biocontainers/quast:5.0.2--py36pl5321hcac48a8_7"
	label "low_computation"

	input:
	path fasta

	output:
	path "*"

	script:
	"""
	quast.py -t $task.cpus -o output *.fasta
	find output -type f -exec mv {} . \\;
	"""
}
