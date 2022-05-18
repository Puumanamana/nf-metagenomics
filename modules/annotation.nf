/*
 * Workflow
 */

workflow ANNOTATE {
	take:
	bins

	main:
	proteins = PRODIGAL(bins).faa

	if (params.amr_detection) {
		pfam_db = Channel.fromPath("$params.pfam_db_dir/*").collect()
		HMMER_HMMSCAN(proteins, pfam_db)
	}

	if (!params.skip_gtdbtk) {
		gtdbtk_db = file(params.gtdbtk_db).isDirectory() ?
			file("$params.gtdbtk_db/db") :
			GTDBTK_DOWNLOAD()

		GTDBTK(
			bins.map{it[1]}.buffer(size: 1000, remainder: true),
			gtdbtk_db
		)
	}
}

/*
 * Processes
 */

process PROKKA {
	tag "$sample"
	publishDir "${params.outdir}/prokka", mode: 'copy'
	container "quay.io/biocontainers/prokka:1.14.6--pl5262hdfd78af_3"
	label "high_computation"

	input:
	tuple val(sample), path(fasta)

	output:
	tuple val(sample), path('prokka_out/*')

	script:
	"""
	prokka $fasta \\
		--outdir prokka_out \\
		--metagenome \\
		--centre X --compliant \\
		--cpus $task.cpus
	"""
}

process PRODIGAL {
	tag "$sample"
	publishDir "${params.outdir}/prodigal", mode: 'copy'
	container "quay.io/biocontainers/prodigal:2.6.3--hec16e2b_4"
	label "medium_computation"

	input:
	tuple val(sample), path(fasta)

	output:
	tuple val(sample), path('*.gbk'), emit: gbk
	tuple val(sample), path('*.faa'), emit: faa

	script:
	"""
	cat $fasta |
	  prodigal -a "${sample}.faa" -p meta \\
	> "genes_${sample}.gbk"
	"""
}

process HMMER_HMMSCAN {
	tag "$sample"
	publishDir "${params.outdir}/resfam", mode: 'copy'
	container "quay.io/biocontainers/hmmer:3.3.2--h1b792b2_1"
	label "high_computation"

	input:
	tuple val(sample), path(fasta)
	path hmm_db

	output:
	tuple val(sample), path("*.tsv"), emit: tsv
	tuple val(sample), path("*.log"), emit: log

	script:
	"""
	db_file=`ls *.hmm`
	hmmscan --tblout "${sample}_resfam.tsv" \$db_file $fasta \
		> "hmmscan_${sample}.log"
	"""
}

process GTDBTK_DOWNLOAD {
	publishDir "$params.gtdbtk_db", mode: "copy"
	container "quay.io/biocontainers/gtdbtk:2.0.0--pyhdfd78af_1"
	label "medium_computation"

	output:
	path "db"

	script:
	"""
	export GTDBTK_DATA_PATH=./db
	mkdir db && download-db.sh
	"""
}

process GTDBTK {
	tag "batch"
	publishDir "${params.outdir}/gtdbtk", mode: 'copy'
	container "quay.io/biocontainers/gtdbtk:2.0.0--pyhdfd78af_1"
	label "high_computation"

	input:
	path fasta
	path db

	output:
	path "taxonomy_gtdbtk/*", emit: taxonomy

	script:
	ext = fasta[0].getExtension()
	"""
	export GTDBTK_DATA_PATH=$db
	gtdbtk classify_wf --genome_dir . --out_dir taxonomy_gtdbtk --extension $ext
	"""
}
