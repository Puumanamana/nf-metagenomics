nextflow.enable.dsl=2

// Functions
include { helpMessage ; saveParams } from "./util.nf"

include { TRIMMING } from './modules/trimming.nf' \
	addParams(outdir: "${params.outdir}/trimming")
include { ASSEMBLY } from './modules/assembly.nf' \
	addParams(outdir: "${params.outdir}/assembly")
include { MERGE_ASSEMBLIES } from './modules/assembly_merging.nf' \
	addParams(outdir: "${params.outdir}/assembly")
include { VIRSORTER2_FILTER } from './modules/virsorter2.nf' \
	addParams(outdir: "${params.outdir}/virsorter2")
include { COVERAGE } from './modules/coverage.nf' \
	addParams(outdir: "${params.outdir}/coverage")
include { BINNING } from './modules/binning.nf' \
	addParams(outdir: "${params.outdir}/binning")
include { ANNOTATE } from './modules/annotation.nf' \
	addParams(outdir: "${params.outdir}/annotation")

workflow WGS {
    take:
	reads
	
    main:
	// Preprocessing
	trimmed = TRIMMING(reads)
	all_reads = trimmed.paired.mix(trimmed.unpaired)

	// Assembly
	contigs = ASSEMBLY(all_reads)
	merged = params.coassembly ? contigs.filt : MERGE_ASSEMBLIES(contigs.filt)

	// Extraction of viral fraction
	if (params.viral) {
		merged = VIRSORTER2_FILTER(merged)
	}

	// Coverage
	cov = COVERAGE(trimmed.paired, merged)

	// Binning
	if (params.binning) {
		BINNING(merged, cov.bam, cov.bai)
		bins = BINNING.out.fasta.flatten()
			.map{[it.getBaseName(), it]}
	} else {
		bins = merged.map{it[1]}.splitFasta(by: 1, file: true)
			.map{[ (it.getBaseName() =~ /([0-9]+)$/ )[0][1], it ]}
	}

	// Annotation
	ANNOTATE(bins)
}

workflow {
	reads = Channel.fromFilePairs(params.reads)
	saveParams()
	WGS( reads )
}
