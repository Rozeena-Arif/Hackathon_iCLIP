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

process GFF_TO_GTF {
	tag  "convert gff to gtf format"
	
	conda 'envs/agat.yml'

        publishDir (
                path: "${baseDir}/assets/annotation/gtf",
                mode: 'copy',
                overwrite: 'true'
                )
        input:
                path(gffFile)
        
        output:
                path 'annotation.gtf', emit: gtf

        script:
                """
                agat_convert_sp_gff2gtf.pl \
                --gff ${gffFile} \
                --gtf_version 3 \
                --output annotation.gtf
                """
}
