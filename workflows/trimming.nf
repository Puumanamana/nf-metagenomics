include { FASTP; TRIMMOMATIC } from '../modules/trimming.nf'
include { READS_QC as QC_BEFORE } from './reads_qc.nf' \
	addParams(outdir: "${params.outdir}/qc_before")
include { READS_QC as QC_AFTER } from './reads_qc.nf' \
	addParams(outdir: "${params.outdir}/qc_after")

