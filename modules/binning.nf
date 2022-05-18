/*
 * Workflow
 */

workflow BINNING {
	take:
	assembly
	bam
	bai

	main:

	cov_summary = JGI_SUMMARIZE_BAM_CONTIG_DEPTHS(
		bam.collect{it[1]},
		bai.collect{it[1]}
	)
	METABAT2(
		assembly.map{it[1]},
		cov_summary
	)

	emit:
	fasta = METABAT2.out.fasta
	tsv = METABAT2.out.tsv
	cov = cov_summary
}

/*
 * Processes
 */

process JGI_SUMMARIZE_BAM_CONTIG_DEPTHS {
	tag "all"
	publishDir "${params.outdir}/metabat2", mode: 'copy'
	container "quay.io/biocontainers/metabat2:2.15--h986a166_1"
	label "low_computation"

	input:
	path bam
	path bai

	output:
	path "depth.txt"

	script:
	"""
	jgi_summarize_bam_contig_depths --outputDepth depth.txt *.bam
	"""
}

process METABAT2 {
	tag "all"
	publishDir "${params.outdir}/metabat2", mode: 'copy'
	container "quay.io/biocontainers/metabat2:2.15--h986a166_1"
	label "medium_computation"

	input:
	path assembly
	path coverage

	output:
	path "assignments.tsv", emit: tsv
	path "bins/*.fa", emit: fasta

	script:
	// s: Minimum bin size (bp)
	// saveCls: Save cluster memberships
	// unbinned: Generate [outFile].unbinned.fa file for unbinned contigs
	"""
	metabat2 -t $task.cpus -i $assembly -a $coverage -o bins/metabat2 \\
	  -s 1 --saveCls --unbinned

	mv bins/metabat2 assignments.tsv
	mv bins/metabat2.{lowDepth,tooShort,unbinned}.fa .
	cat metabat2.*.fa | awk '/^>/ {fname=substr(\$1, 2)".fa"} {print > "bins/"fname}'
	"""
}
