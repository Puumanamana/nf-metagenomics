include { JGI_SUMMARIZE_BAM_CONTIG_DEPTHS; METABAT2 } from '../modules/binning.nf'

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
