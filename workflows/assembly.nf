include { MEGAHIT; SPADES; QUAST } from '../modules/assembly.nf'
include { SEQTK_SEQ as MIN_LENGTH_FILTER } from '../modules/seqtk.nf'
include { SEQTK_SUBSEQ } from '../modules/seqtk.nf'
include { CONCAT_ASSEMBLIES; PAIRWISE_BLAST; ANI_CLUSTERING } from '../modules/assembly_merging.nf'

