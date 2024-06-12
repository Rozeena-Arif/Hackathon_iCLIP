include {CONCAT} from './modules/concat_genome.nf'
include {GFF_TO_GTF} from './modules/concat_genome.nf'

workflow
{
	CONCAT(fastaDir=params.refGenomeDir)
	GFF_TO_GTF(gffFile="${params.annotation}${params.gff_file}")
}

//
