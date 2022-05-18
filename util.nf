// Help function for metagenomics pipeline

def helpMessage() {
	log.info"""
	===================================
	WGS-pipeline
	===================================
	Usage:

	---------------------------------- Mandatory arguments ----------------------------------------
	--reads         Path to input data (glob pattern)

	---------------------------------- Optional arguments ----------------------------------------
	--outdir        Path to output directory. Default: "./outputs-wgs"

	[Trimming]
	--trimmer            Trimming tool to use. Choices: [FASTP] or Trimmomatic
	--adapters           Path to adapter file used for sequencing.
						 If empty, auto-detection by Fastp
	--cut_mean_quality   Trim read if sequence quality drops below this value in a sliding window
						 from 3' to 5'. Default: [25]
	--cut_window_size    Window size for quality trimming. Default: [4]
	--trim_front1        FASTP: Trim bases at the 5' end of forward read.
						 Trimmomatic: Trim bases at the 5' end of both reads. Default: 0
	--trim_front2        Trim bases at the 5' end of reverse read (FASTP only). Default: 0
	--trim_tail1         FASTP: Trim bases at the 3' end of forward read.
						 Trimmomatic: Trim bases at the 3' end of both reads. Default: 0
	--trim_tail2         Trim bases at the 3' end of reverse read (FASTP only). Default: 0

	[Assembly]
	--assembler          Any spades assembler, or megahit. Use megahit for large datasets.
						 Default: metaspades
	--coassembly         Pool samples before assembly. Default: False

	[Contig filtering]
	--min_contig_len     Default: 2000 (bp)
	--viral              Keep only viral contigs with virsorter2

	[Contig coverage]
	--min_mapq           Min read mapping quality. Default: 30

	[Annotation]
	--amr_detection      Look for AMR genes. Default: False
	--pfam_db_dir        HMM database for AMR gene detection. Default is the resfam database
						 Available at http://dantaslab.wustl.edu/resfams/Resfams-full.hmm.gz
	""".stripIndent()
}

// Show help message
params.help = false
if (params.help){
	helpMessage()
	exit 0
}

def saveParams() {
	file(params.outdir).mkdir()
	File f = new File("${params.outdir}/parameters_summary.log")
	f.write("====== Parameter summary =====\n")
	params.each { key, val -> f.append("\n${key.padRight(50)}: ${val}") }
}
