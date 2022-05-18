include { SEQTK_SUBSEQ } from './seqtk.nf'

/*
 * Workflow
 */

workflow MERGE_ASSEMBLIES {
	take: assemblies

	main:
	// concatenate all fasta files
	concat = CONCAT_ASSEMBLIES(assemblies.collect{it[1]})
	// pairwise comparison of all contigs
	blast = PAIRWISE_BLAST(concat)
	// ANI-based clustering
	clusters = ANI_CLUSTERING(concat.join(blast)).clusters
	// Keep representative (=longest) contig of each cluster
	merged = SEQTK_SUBSEQ(concat.join(clusters))
		.mix(concat).first() // if no clusters, keep concat

	emit:
	merged=merged
}

/*
 * Processes
 */

process CONCAT_ASSEMBLIES {
	tag "all"
	label "low_computation"
	label "py_script"

	input:
	path assemblies

	output:
	tuple val("pooled"), path("concat.fasta")

	script:
	def fasta = assemblies.join("','")
	"""
	#!/usr/bin/env python

	from pathlib import Path
	from Bio.SeqIO.FastaIO import SimpleFastaParser

	with open("concat.fasta", "w") as writer:
		for fa in ['$fasta']:
			with open(fa) as reader:
				for (title, seq) in SimpleFastaParser(reader):
					writer.write(f">{title}.{Path(fa).stem}\\n{seq}\\n")
	"""
}

process PAIRWISE_BLAST {
	tag "$sample"
	container "quay.io/biocontainers/blast:2.12.0--hf3cf87c_4"
	label "high_computation"

	input:
	tuple val(sample), path(assembly)

	output:
	path "blast.tsv", optional: true

	script:
	"""
	makeblastdb -in $assembly -dbtype nucl -out db
	blastn -query $assembly -db db -outfmt '6 std qlen slen' -max_target_seqs 10000 -perc_identity 90 -out - -num_threads $task.cpus |
		awk '\$1 != \$2' > blast.tsv
	[ ! -s blast.tsv ] && rm -f blast.tsv || echo "Hits found"
	"""
}

process ANI_CLUSTERING {
	tag "$sample"
	label "py_script"
	label "medium_computation"

	input:
	tuple val(sample), path(assembly), path(blast)

	output:
	tuple val(sample), path("clusters*.tsv"), emit: clusters
	tuple val(sample), path("ani*.tsv"), emit: ani

	script:
	"""
	anicalc.py -i $blast -o ani_${sample}.tsv
	aniclust.py --fna $assembly --ani ani_${sample}.tsv --out clusters_${sample}.tsv \\
	  --min_ani 95 --min_tcov 85 --min_qcov 0
	"""
}
