include { SEQTK_SUBSEQ } from '../modules/seqtk.nf'
include { VIRSORTER2_SETUP; VIRSORTER2_RUN } from '../modules/virsorter2.nf'

workflow VIRSORTER2_FILTER{
	take: contigs

	main:
	db = VIRSORTER2_SETUP()
	vids = VIRSORTER2_RUN(contigs, db).ids

	filtered = SEQTK_SUBSEQ(contigs.join(vids))

	emit:
	filtered = filtered
}
