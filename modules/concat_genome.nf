process CONCAT {
	tag "concating the fasta files"

	publishDir (
		path: "${baseDir}/assets/reference_genome",
		mode: 'copy',
		overwrite: 'true'
		)
	input:
		path(fastaDir)

	output: 
		path 'all_sequences.fa', emit: ref_genome

script:
	"""
	cat ${fastaDir}/*.fa > all_sequences.fa
	"""
}

