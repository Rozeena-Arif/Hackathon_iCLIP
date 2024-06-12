include {CONCAT} from './modules/concat_genome.nf'
include {GFF_TO_GTF} from './modules/gff_to_gtf.nf'

workflow
{
	CONCAT(fastaDir=params.refGenomeDir)
	GFF_TO_GTF(gffFile=params.annotation_dir/params.gff_file)
}
