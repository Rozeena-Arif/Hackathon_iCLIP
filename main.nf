include {CONCAT} from './modules/concat_genome.nf'

workflow
{
	CONCAT(fastaDir=params.refGenomeDir)
}
