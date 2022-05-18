include { SEQTK_SUBSEQ } from './seqtk.nf'

/*
 * Workflow
 */

workflow VIRSORTER2_FILTER{
	take: contigs

	main:
	db = VIRSORTER2_SETUP()
	vids = VIRSORTER2_RUN(contigs, db).ids

	filtered = SEQTK_SUBSEQ(contigs.join(vids))

	emit:
	filtered = filtered
}

/*
 * Processes
 */

process VIRSORTER2_SETUP {
	tag "db"
	label "medium_computation"
	container "docker://jiarong/virsorter:2.2.3"

	output:
	path "virsorter2_db", emit: db

	script:
	"""
	virsorter setup -d virsorter2_db -j $task.cpus
	"""
}

process VIRSORTER2_RUN {
	tag "$sample"
	label "medium_computation"
	container "docker://jiarong/virsorter:2.2.3"

	input:
	tuple val(sample), path(fasta)
	path virsorter2_db

	output:
	tuple val(sample), path("virsorter2-*.txt"), optional: true, emit: ids
	tuple val(sample), path("vs2_out/final-viral-score.tsv"), emit: scores

	script:
	ids_file = "virsorter2-${sample}.txt"
	"""
	virsorter run -w vs2_out -i $fasta -j $task.cpus -d $virsorter2_db
	tail -n+2 final-viral-score.tsv | cut -f1 | cut -d'|' -f1 > $ids_file

	[ ! -s $ids_file ] && rm -f $ids_file || echo "viruses found"
	"""
}
