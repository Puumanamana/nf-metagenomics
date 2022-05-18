process SEQTK_SEQ {
	tag "$sample"
	publishDir "${params.outdir}/seqtk", mode: 'copy'
	container "quay.io/biocontainers/seqtk:1.3--h7132678_4"
	label "low_computation"

	input:
	tuple val(sample), path(fasta)

	output:
	tuple val(sample), path('*.fasta'), emit: fasta

	script:
	"""
	seqtk seq -L $params.min_contig_len $fasta \\
	   > ${sample}_gt${params.min_contig_len}.fasta
	"""
}

process SEQTK_SUBSEQ {
	tag "$sample"
	publishDir "${params.outdir}/seqtk", mode: 'copy'
	container "quay.io/biocontainers/seqtk:1.3--h7132678_4"
	label "low_computation"

	input:
	tuple val(sample), path(fasta), path(list)

	output:
	tuple val(sample), path("*.fasta")

	script:
	"""
	seqtk subseq $fasta $list > ${fasta.getSimpleName()}_subset.fasta
	"""
}
