##########################################
##########################################
############ INPUT FUNCTIONS #############
##########################################
##########################################

get_gr_and_bs <- function(filepath, protein, condition){
  
  bed <- rtracklayer::import.bed(filepath)
  bed <- granges(bed)
  seqlevelsStyle(bed) <- "UCSC"
  
  bed$id <- paste(seqnames(bed), ranges(bed), as.character(strand(bed)), sep='_' )
  
  bed$condition <- condition
  bed$protein <- protein
  
  reduced <- bed %>% GenomicRanges::reduce()
  reduced$id <- paste(seqnames(reduced), ranges(reduced), as.character(strand(reduced)), sep='_' )
  binding_sites <- get_binding_sites(genome_annotations, reduced) %>% unique()
  intron_binding_sites <- intron_binding(data = reduced, intron_annotations = introns)
  all_binding_sites <- merge(binding_sites, intron_binding_sites, all.x = TRUE, all.y = TRUE) %>% fill(10:12) %>% 
    add_column(condition = condition, protein = protein)
  
  
  list(gr = bed, bs = all_binding_sites)
}

get_bs_count_info <- function(bs,gr,var=NA, background =NA){
  protein <- bs$protein %>% unique()
  condition <- bs$condition %>% unique()
  
  n_target_genes <- bs %>% as_tibble() %>% pull(gene_name) %>% unique() %>% length()
  n_unique_binding_sites <- gr %>% length()
  n_unique_binding_sites_nonoverlapping <- gr %>% GenomicRanges::reduce() %>% length()
  
  count_tbl <- tibble(type = factor(c("n_unique_binding_sites", 
                                      "n_unique_binding_sites_nonoverlapping", 
                                      "n_target_genes"), 
                                    levels = c("n_unique_binding_sites", 
                                               "n_unique_binding_sites_nonoverlapping", 
                                               "n_target_genes")), 
                      count = c(n_unique_binding_sites, n_unique_binding_sites_nonoverlapping, n_target_genes),
                      protein = protein,
                      condition = condition)
  
  if(!is.na(var)){
    count_tbl$var <- var
  }
  
  if(!is.na(background)){
    count_tbl$background <- background
  }
  
  count_tbl
  
} 

####### FOREGROUND INPUT FUNCTIONS ####### 
get_iclip_input_files <- function(iclip_bed_filepath, iclip_bw_plus_filepath, iclip_bw_minus_filepath){
  
  #read in bed file 
  #this is not actually a bed file in this case but is treated as one
  #it contains the genome coordinates for iCLIP hits
  bed <- read.delim(iclip_bed_filepath, header = FALSE)
  
  #select relevant columns from input 'bed' to created a granges object
  bed <- bed[-1,] %>% 
    dplyr::select(V1, V2, V3, V6) %>% 
    dplyr::rename(seqnames = V1, start = V2, end = V3, strand = V6) %>% 
    as_granges()
  
  #convert the chromosome name format to be 'chr1' instead of '1'
  #this is important for future steps/ functions which use the UCSC style
  seqlevelsStyle(bed) <- "UCSC"
  
  #combined any ranges that might be overlapping
  bed_reduced <- bed %>% GenomicRanges::reduce()
  
  #create ids for each unique iCLIP hit 
  #these contain the chromosome name, the start and end coordinates, and the strand 
  ids <- paste(seqnames(bed_reduced), ranges(bed_reduced), as.character(strand(bed_reduced)), sep='_' )
  #add these ids to the granges object
  bed_reduced$id <- ids
  
  #import bigwig files, which have information on level of signal at each position 
  #there are separate bigwig files for the plus and minus strands so they are imported separately
  bw_p <- import.bw(iclip_bw_plus_filepath)
  seqlevelsStyle(bw_p) <- "UCSC"
  
  bw_m <- import.bw(iclip_bw_minus_filepath)
  seqlevelsStyle(bw_m) <- "UCSC"
  
  iclip_files <- list(bed = bed_reduced, bigwig_plus = bw_p, bigwig_minus = bw_m)
  
  return(iclip_files)
  
}

bw_bed_overlap <- function(bed, bw){  
  bed_gr_bwp_overlaps <- findOverlaps(bw, bed)
  
  bigwig_overlap_summary <- data.frame(chromosome = seqnames(bw[queryHits(bed_gr_bwp_overlaps)]),
                                       range = ranges(bw[queryHits(bed_gr_bwp_overlaps)]),
                                       id = bed[subjectHits(bed_gr_bwp_overlaps)]$id,
                                       score = score(bw[queryHits(bed_gr_bwp_overlaps)])) %>% 
    distinct()
  
  
  bigwig_normalised <- bigwig_overlap_summary %>% 
    #this takes the start and end of each range and stretches it out
    dplyr::mutate(position = map2(range.start, range.end, ~ seq(from = .x, to = .y))) %>% 
    #this selects only the columns we're interested in
    unnest(c(id, position, score))  %>%  
    #to normalise the data, we divide them all by the sum of scores and multiply by 1,000,000 to make the numbers more user-friendly
    dplyr::mutate(normalised = score) 
  
  return(bigwig_normalised)
}

iCLIP_prediction_foreground_peaks <- 
  function(bed, bigwig_plus, bigwig_minus){
    
    #the bigwig signal information needs to be matched up with the bed coordinates
    #this is achieved using the function bw_bed_overlap
    #it is done separately for plus and minus strands
    bw_bed_p <- bw_bed_overlap(bed, bigwig_plus)
    bw_bed_m <- bw_bed_overlap(bed, bigwig_minus)
    
    #the strand is added as a column to each of the tibbles generated
    bw_bed_p$strand <- "+"
    bw_bed_m$strand <- "-"
    
    #the dataframes are combined into a single tibbles
    bw_bed <- rbind(bw_bed_p, bw_bed_m)
    
    #the peak of each binding site is identified by identifying the position with the highest signal
    bw_bed_peaks <- bw_bed %>% 
      group_by(id) %>% 
      dplyr::mutate(peak = max(normalised)) %>% 
      dplyr::filter(peak == normalised) 
    
    #a new granges object is created where the ranges are only the peak position
    #this can be used as the starting point for creating inputs for sequence prediction
    peak_granges <- GRanges(seqnames = bw_bed_peaks$chromosome, 
                            ranges = IRanges(start = bw_bed_peaks$position, end = bw_bed_peaks$position), 
                            strand = bw_bed_peaks$strand,
                            id = bw_bed_peaks$id,
                            peak_position = bw_bed_peaks$position,
                            normalised_score = bw_bed_peaks$normalised)
    
    return(peak_granges)
  }

iCLIP_prediction_foreground_fasta <- 
  function(peak_granges, width, genome, output = "fasta"){
    
    
    #extend the coordinates of each hit to be the desired width for motif prediction
    gr_for_fasta <- peak_granges %>% flank(width/2, both = TRUE)
    
    if(output == "fasta"){
    #some of the iCLIP datasets contain binding sites within viral genomes
    #these can't be extracted from the human genome so should be removed in advance
    #they can be manually added later if desired for motif prediction
    #first, chromosome names are extracted from the genome object
    genome_names <- seqnames(genome)
    #these are used to filter the peaks 
    filtered_peaks <- gr_for_fasta %>% 
      plyranges::filter(seqnames %in% genome_names)
    #the viral genome will still be present in the levels of the the chromosome column
    #it needs to be removed for sequence extraction to work
    filtered_peaks  <- keepSeqlevels(filtered_peaks , unique(seqnames(filtered_peaks)@values) )
    
    fasta <- Biostrings::extractAt(genome, filtered_peaks)
    
    names(fasta) <- filtered_peaks$id
    return(fasta)
    
    }else{return(gr_for_fasta)}
    
    
  }

####### BINDING SITE FUNCTIONS ####### 

#this function establishes where the binding sites identified in the eclip experiments land in the genome
#it takes the genome annotations and the imported eclip data as inputs and returns a dataframe with binding site information
get_binding_sites <- function(annotations, data){
  #this uses findOverlaps function from GRanges, which looks for overlaps in the coordinates in 2 GRanges objects
  #it returns a table with corresponding row numbers of overlapping regions 
  #queryHits are the row numbers from the eclip dataset and subjectHits and the row numbers from the genome annotation
  overlaps <- findOverlaps(data, annotations)
  #the row numbers supplied by findOverlaps can be used to extract information from Granges 
  #this information can be compiled into a dataframe, which is then returned 
  overlap_summary <- data.frame(query = queryHits(overlaps),
                                chromosome = seqnames(data[queryHits(overlaps)]),
                                range = ranges(data[queryHits(overlaps)]),
                                strand = strand(data[queryHits(overlaps)]),
                                id = data[queryHits(overlaps)]$id,
                                gene_id = annotations[subjectHits(overlaps)]
                                $gene_id,
                                transcript_id = annotations[subjectHits(overlaps)]
                                $transcript_id,
                                gene_name = annotations[subjectHits(overlaps)]
                                $gene_name,
                                gene_biotype = annotations[subjectHits(overlaps)]
                                $gene_biotype,
                                type = annotations[subjectHits(overlaps)]$type,
                                exon_number = annotations[subjectHits(overlaps)]$exon_number,
                                exon_range = ranges(annotations[subjectHits(overlaps)])) 
  overlap_summary <- distinct(overlap_summary)
  return(overlap_summary)
}


#this function takes the eclip data and annotation file generated for introns and returns a dataframe with all the binding sites present in introns
intron_binding <- function(data, intron_annotations){
  #as above, this creates a table with references to the rows in which there are matches between the query data and the annotations
  overlaps <- findOverlaps(data, intron_annotations)
  #this information can be used to generate a dataframe     
  intron_overlaps <- data.frame(query = queryHits(overlaps),
                                chromosome = seqnames(data[queryHits(overlaps)]),
                                range = ranges(data[queryHits(overlaps)]),
                                strand = strand(data[queryHits(overlaps)]),
                                id = data[queryHits(overlaps)]$id,
                                gene_id = intron_annotations[subjectHits(overlaps)]
                                $gene_id,
                                type = "intron")  
  return(intron_overlaps)
}

iCLIP_prediction_background_fasta <- 
  function(peak_granges, width, genome_annotations, intron_annotations, genome){
    
    genome_names <- seqnames(genome)
    filtered_peaks <- peak_granges %>% 
      plyranges::filter(seqnames %in% genome_names)
    filtered_peaks  <- keepSeqlevels(filtered_peaks , unique(seqnames(filtered_peaks)@values) )
    
    #to generate gene and region-matched input, it is first necessary to determine the genes/regions that binding sites are occurring in
    #first, identify binding sites in all regions other than introns using normal genome annotation
    binding_sites <- get_binding_sites(genome_annotations, filtered_peaks) %>% unique()
    #next, identify binding sites in introns using custom intron annotation
    #they are done separately because introns are not annotated in the genome annotation object
    intron_binding_sites <- intron_binding(data = keepStandardChromosomes(filtered_peaks), 
                                           intron_annotations = intron_annotations)
    #combine these two and fill any gaps in intron annotation with information from genome annotation
    all_binding_sites <- merge(binding_sites, intron_binding_sites, all.x = TRUE, all.y = TRUE) %>% fill(10:12)
    
    #we start by generating the background for sequences that are non-coding 
    #this is done because there are some sequences that will end up being annotated both as intron and as a non-coding gene
    #to ensure background for these hits is not generated twice, we only consider them as non-coding genes
    counts <- all_binding_sites %>% 
      dplyr::filter(gene_biotype != "protein_coding") %>% 
      distinct(id, .keep_all = T) %>%
      count(gene_id) %>% 
      add_column(type = "exon") 
    
    #filter the genome annotation to only include genes for which background is required
    #there might be multiple transcripts for the same gene
    anno <- as_tibble(genome_annotations) %>% 
      dplyr::filter(gene_id %in% counts$gene_id & type == "exon") %>%
      dplyr::select(seqnames, strand, start, end, gene_id, transcript_id, exon_number) %>% as_granges()
    
    #extract the sequence for these genes
    sequences <- Biostrings::extractAt(genome, anno)
    
    #create a table with the gene and transcript ids, exon number and actual sequence
    exon_sequences <- tibble(
      gene_id     = anno$gene_id,
      tx_id       = anno$transcript_id,
      exon_number = as.numeric(anno$exon_number),
      sequence    = as.character(sequences)
    )
    
    #exon number is used to arrange the exons for each transcript
    #this means that transcripts can be pasted together
    transcript_sequences <- exon_sequences %>% 
      group_by(tx_id) %>% 
      arrange(exon_number) %>% 
      dplyr::summarise(sequence = str_c(sequence, collapse = "")) %>% 
      as_tibble()
    
    #
    transcript_sequences <- exon_sequences %>% 
      dplyr::select(gene_id, tx_id) %>% 
      distinct() %>% 
      left_join(transcript_sequences, by = "tx_id") %>% 
      dplyr::select(-tx_id) %>% 
      add_column(type = "gene")
    
    transcript_sequences <-  transcript_sequences %>% 
      dplyr::filter(width(sequence) > width) 
    
    counts <- counts %>% dplyr::filter(gene_id %in% transcript_sequences$gene_id)
    
    nc_background <- map2_df(counts$gene_id, counts$n, function(gene, freq){
      
      seqs <- transcript_sequences  %>% 
        dplyr::filter(gene_id == gene) %>% 
        pull(sequence)
      
      background <- sample_sequences(seqs, freq, width)
      tibble(gene = gene, type = "exon", sequence = background)
    }
    ) 
    
    ids_nc <- all_binding_sites %>% 
      dplyr::filter(gene_id %in% nc_background$gene_id) %>% 
      distinct(id, .keep_all = T) 
    
    remaining <- all_binding_sites %>% dplyr::filter(!(id %in% ids_nc$id))
    remaining_noint <- remaining %>% dplyr::filter(type != "intron")
    remaining_int <- remaining %>% dplyr::filter(type == "intron") 
    
    # Get number of counts per each gene and region
    counts <- remaining_noint %>% 
      filter(type == "CDS" | type == "three_prime_utr" | type == "five_prime_utr") %>% 
      distinct(id, .keep_all = TRUE) %>% 
      group_by(gene_id) %>% count(type) %>% 
      ungroup() 
    
    genes <- counts$gene_id
    types <- counts$type
    
    anno <- as_tibble(genome_annotations) %>% 
      dplyr::filter(gene_id %in% genes, type %in% types) %>%
      dplyr::select(seqnames, strand, start, end, type, gene_id, transcript_id, exon_number) %>% as_granges()
    
    # Extract sequences
    sequences <- Biostrings::extractAt(genome, anno)
    
    # Build table with sequence per exon
    exon_sequences <- tibble(
      gene_id     = anno$gene_id,
      tx_id       = anno$transcript_id,
      type        = anno$type,
      exon_number = as.numeric(anno$exon_number),
      sequence    = as.character(sequences)
    )
    
    # Build transcript units from exons
    transcript_sequences <- exon_sequences %>% 
      group_by(tx_id, type) %>% 
      arrange(exon_number) %>% 
      dplyr::summarise(
        sequence = str_c(sequence, collapse = ""), 
        .groups  = "drop"
      )
    
    # Attach gene id to transcript table
    transcript_sequences <- exon_sequences %>% 
      dplyr::select(gene_id, tx_id) %>% 
      distinct() %>% 
      left_join(transcript_sequences, by = "tx_id") %>% 
      dplyr::select(-tx_id)
    
    protein_coding_bkg <- pmap_df(list(counts$gene_id, counts$type, counts$n), 
                                  sample_transcripts,
                                  transcript_sequences, 
                                  width = width)
    
    
    counts <- as_tibble(remaining_int) %>% 
      distinct(id, .keep_all = TRUE) %>% 
      count(gene_id) 
    counts$type <- "intron"
    
    
    anno <- as_tibble(introns) %>% 
      dplyr::filter(gene_id %in% counts$gene_id) %>%
      dplyr::select(seqnames, strand, start, end, gene_id) %>% as_granges()
    
    sequences <- Biostrings::extractAt(genome, anno)
    
    
    seq_table <- tibble(
      gene_id = anno$gene_id,
      sequence = as.character(sequences),
      type = "intron"
    )
    
    seq_table <- seq_table %>% dplyr::filter(width(sequence) > width)
    
    counts <- counts %>% dplyr::filter(gene_id %in% seq_table$gene_id)
    
    
    intron_background <- pmap_df(list(counts$gene_id, counts$type, counts$n), 
                                 sample_transcripts,
                                 seq_table, 
                                 width = width) 
    
    
    
    all_bkg <- rbind(intron_background, protein_coding_bkg, nc_background)
    
    background_seqs <- DNAStringSet(all_bkg$sequence)
    names(background_seqs) <- all_bkg$gene %>% make.unique()
    
    return(background_seqs)
    
  }


get_random_subseq <- function(sequence, width) {
  pos <- sample(1:(nchar(sequence) - width), 1)
  seq <- str_sub(sequence, pos, pos + width - 1)
  seq
}

sample_sequences <- function(sequences, n, width) {
  seqs  <- sequences[sample(1:length(sequences), n, replace=TRUE)]
  bg_seqs <- map_chr(seqs, get_random_subseq, width)
  bg_seqs
}

sample_transcripts <- function(gene, type, n, tx_table, width = 50) {
  
  # Get sequences corresponding to gene and type form tx_table
  sequences <- tx_table %>% 
    dplyr::filter(gene_id == gene, type == type) %>% 
    pull(sequence)
  
  # Exclude sequences that are too short (shorter than requested with)
  sequences <- sequences[nchar(sequences) >= width]
  if(is_empty(sequences)) {
    warning(str_c("No long enough sequence available for gene '", gene, 
                  "' type '", type, "' (width requested is ", width, ")"))
    return(NA)
  }
  
  # Sample
  background <- sample_sequences(sequences, n, width)
  tibble(gene = gene, type = type, sequence = background)
}

########## BACKGROUND FUNCTIONS ##########



##########################################
##########################################
####### MOTIF ANALYSIS FUNCTIONS #########
##########################################
##########################################

motif_coverage <- function(sequences, motifs, threshold = 70){
  
  PURA_HIV_scan_red <- scan_sequences(motifs, sequences, threshold = 0.01) %>%
    as_tibble() %>% 
    dplyr::filter(score.pct > threshold)
  
  
  motif_seqs_sep <- PURA_HIV_scan_red %>% 
    group_split(motif.i) %>% 
    map(~ pull(.x, sequence)) %>% 
    map(~ unique(.x))
  
  PURA_HIV_scan_red_co <- PURA_HIV_scan_red %>% 
    dplyr::mutate(co_occurrence = case_when(
      sequence %in% motif_seqs_sep[[1]] & sequence %in% motif_seqs_sep[[2]] ~ "both",
      sequence %in% motif_seqs_sep[[1]] & !(sequence %in% motif_seqs_sep[[2]]) ~ "motif_1",
      sequence %in% motif_seqs_sep[[2]] & !(sequence %in% motif_seqs_sep[[1]]) ~ "motif_2"
    ))
  
  all_ids <- names(sequences)
  
  none <- setdiff(all_ids, PURA_HIV_scan_red_co$sequence)
  
  counts <- count(PURA_HIV_scan_red_co$co_occurrence) %>% 
    add_row(x = "none", freq = length(none)) %>% 
    dplyr::mutate(proportion = freq/sum(freq)) %>% 
    dplyr::rename(motif_occurrence = x)
  
  return(counts)
  
}

motif_scan_with_cooccurrence <- function(sequence, motifs, threshold){
  
  scan <- scan_sequences(motifs, sequence, 0.001) %>% 
    as_tibble() %>% dplyr::filter(score.pct >= threshold)
  
  scan_split <- scan %>% 
    group_split(motif.i) %>% 
    map(~ pull(.x, sequence)) %>% 
    map(~ unique(.x))
  
  scan_with_co <- scan %>% 
    dplyr::mutate(co_occurrence = case_when(
      sequence %in% scan_split[[1]] & sequence %in% scan_split[[2]] ~ "both",
      sequence %in% scan_split[[1]] & !(sequence %in% scan_split[[2]]) ~ "motif_1",
      sequence %in% scan_split[[2]] & !(sequence %in% scan_split[[1]]) ~ "motif_2"
    ))
  
  return(scan_with_co)
}

motif_cooccurrence_relative_positions <- function(scan_with_cooccurrence){
  
  new <- scan_with_cooccurrence %>% 
    dplyr::filter(co_occurrence == "both") %>% 
    group_split(motif.i) 
  
  new[[1]] %>% 
    left_join(new[[2]], by = "sequence") %>% 
    dplyr::mutate(difference = start.x - start.y)  
  
}


##########################################
##########################################
########## STRUCTURE FUNCTIONS ###########
##########################################
##########################################


#this generates input sequences of the desired length centred around motif occurrences 
#it takes the motif, input fasta, granges used to generate that fasta,
#score cutoff (for motif occurrences), genome, genome annotation and desired sequence width as input
#it returns a DNAStringSet of the sequences to be folded
input_for_structure_preference <-  function(motif, input_fasta, granges_for_fasta, score_cutoff = 60, genome, genome_annotation, width = 200){
  
  #scan input sequences for motif occurrences and narrow down based on desired threshold
  scan <- scan_sequences(motif, input_fasta) %>% 
    as.tibble() %>% dplyr::filter(score.pct > score_cutoff) %>% 
    dplyr::rename(id = sequence)
  
  #add the granges used to generate the input fasta to the scan output table
  #this relies in merging based on binding site ids
  #it is important so that motif sites can be related back to genome coordinates
  with_motif <- scan %>% 
    left_join(as.tibble(granges_for_fasta), by = "id")
  
  #get length of motif
  motif_length <- motif["consensus"] %>% width()
  
  #put motif coordinates into genome coordinates for plus and minus strands
  
  with_motif_p <- with_motif %>%
    dplyr::filter(strand.y == "+") %>%
    dplyr::mutate(start = start.y + (start.x-1)) %>%
    dplyr::mutate(end = start + motif_length-1) %>%
    dplyr::select(strand.y, seqnames, start, end, id) %>%
    dplyr::rename(strand = strand.y) %>%
    dplyr::filter(seqnames != "HIV-Nef-mCh") %>%
    dplyr::mutate(seqnames = droplevels(seqnames))
  
  
  with_motif_m <- with_motif %>%
    dplyr::filter(strand.y == "-") %>%
    dplyr::mutate(start = end - (stop-1)) %>%
    dplyr::mutate(end = end - (start.x-1)) %>%
    dplyr::select(strand.y, seqnames, start, end, id) %>%
    dplyr::rename(strand = strand.y) %>%
    dplyr::filter(seqnames != "HIV-Nef-mCh") %>%
    dplyr::mutate(seqnames = droplevels(seqnames))
  
  #combine plus and minus tables
  #convert them to granges and extend the width to desired length
  for_secondary_structure <- rbind(with_motif_p, with_motif_m) %>% as_granges() %>% flank(width=width/2, both=TRUE)
  
  #extract sequence
  seqs_for_ss <- Biostrings::extractAt(genome, for_secondary_structure)
  
  #name by ids of the sequence that the motif was identified in 
  names(seqs_for_ss) <- for_secondary_structure$id
  
  return(seqs_for_ss)
}

### helper functions for scrambling ###
explode <- function(dna_seq) {
  substring(text = dna_seq,
            first = seq(1, str_length(dna_seq), 1),
            last  = seq(1, str_length(dna_seq), 1))
}

scramble_dna <- function(dna_seq) {
  dna_seq %>%
    explode() %>%
    sample(x = .,
           size = str_length(dna_seq),
           replace = FALSE) %>%
    str_c(collapse = "")
}

########################################


#this function generates scrambled background 
#it takes the original sequences as input and generates a set number of scrambles for each
#it returns a DNAStringSet of the scrambles
create_scrambled_background <- function(input_sequence, number_of_scrambles = 10){
  
  scramble <- map(as.character(input_sequence), function(x){
    map(seq(1:number_of_scrambles), ~ scramble_dna(x))
  } )
  
  scramble_ss <- scramble %>% unlist() %>%  DNAStringSetList() %>% unlist()
  
  return(scramble_ss)
}

#runs rnafold in input file (needs to already be saved)
#reads and returns output file
run_rnafold <- function(input_filename){
  
  output_filename <- input_filename %>% 
    str_remove(".fa") %>% 
    paste0(., "_folded.fa")
  
  file.create(output_filename)
  
  system(paste0("rnafold --noPS -o", output_filename, " < ", input_filename))
  
  folded_output <- read.table(output_filename, sep = "\t") %>% as_tibble()
  
  return(folded_output)
}

#this function has different options for processing rnafold output
#id_db will return the id and dot bracket for each sequence (this is the required input for forgi)
#db will return just the dotbracket for each sequence
#mfed will return the numeric value for each sequence
process_rnafold_output <- function(rna_fold_output, format = c("id_db", "db", "mfed")){
  
  if(format == "id_db"){
    split(rna_fold_output,rep(1:(nrow(rna_fold_output)/3),each=3)) %>% 
      map2(., seq(from = 1, to = length(.)), function(x, y){
        len <- x[2,] %>% pull() %>% RNAString() %>% length()
        db <- x[3,] %>% strtrim(len) %>% as.tibble() %>% dplyr::rename(V1 = 1)
        id <- ifelse(x[1,] == ">", paste0(">", y), as.character(x[1,])) %>% tibble(V1 = .)
        rbind(id ,db)
      }) %>% ldply(rbind) %>% dplyr::select(V1)
  }else if(format == "db"){
    split(rna_fold_output,rep(1:(nrow(rna_fold_output)/3),each=3)) %>% 
      map(function(x){
        len <- x[2,] %>% pull() %>% RNAString() %>% length()
        x[3,] %>% strtrim(len) %>% as.tibble() %>% dplyr::rename(V1 = 1)
      }) %>% ldply(rbind) %>% dplyr::select(V1)
  }else if(format == "mfed"){
    split(rna_fold_output,rep(1:(nrow(rna_fold_output)/3),each=3)) %>% 
      map(function(x){
        splitmfed <- x[3,] %>% str_split(" ")
        splitmfed[[1]][[2]] %>% 
          str_remove(., "\\(") %>% str_remove(., "\\)") %>% 
          as.numeric() %>% tibble(MFED = .)
      }) %>% ldply(rbind) %>% dplyr::select(MFED)
    
  }
  
}

#this function runs forgi (which assigns structural features to dotbrackets)
#takes as input the filename of saved dotbracket file 
#it returns the forgi output file
run_forgi <- function(input_filename){
  
  directory_filename <- input_filename %>% 
    str_remove(".dotbracket") %>% 
    paste0(., "_forgi_output")
  
  #if the output filename already exists, it will append new output to it, rather than overwriting it
  #to avoid this, this function will remake the file each time
  make.dir2 <- function(fp) {
    if(!file.exists(fp)) {  # If the folder does not exist, create a new one
      #make.dir(dirname(fp))
      dir.create(fp)
    } else {   # If it existed, delete and replace with a new one  
      unlink(fp, recursive = TRUE)
      dir.create(fp)
    }
  } 
  
  make.dir2(directory_filename)
  
  #paste input to be run from command line 
  system_input <- paste0("/Library/Frameworks/Python.framework/Versions/3.8/bin/rnaConvert.py ", 
                         input_filename, " -T element_string --filename ", directory_filename, "/forgi_output")
  
  #run forgi from command line
  system(system_input)
  
  #forgi output generates one file per input sequence 
  #here, we load the names of all these files into R
  filenames_forgi <- list.files(path = directory_filename, full.names = T)
  
  #these names can be used to read all the files in 
  forgi_out <- map(filenames_forgi, ~ read.delim(.x, header = FALSE))
  
  return(forgi_out)
}

#this function takes forgi output (in the form or a list of tables) as input
#it returns a table with a tally of the proportion of each feature at each position in the sequence 
process_forgi_output <- function(forgi_output){
  
  output_length <- forgi_output[[1]][2,] %>% width()
  
  #combine all the outputs into one table
  #separate out the character string so that there is one position per column
  col_per_position <- forgi_output %>% map(~ as.tibble(.x[2,])) %>% 
    ldply(rbind) %>% 
    separate(col = 1, into = as.character(seq(from = 0, to = output_length)), "") 
  
  #get the tally of each feature at each position
  tally_per_position <- col_per_position %>% map(function(x){
    .x <-  paste(x, collapse="")
    total <- width(.x)
    tibble(
      feature = 
        c("unpaired_fiveprime", "stem", "interior_loop", "multiloop_segment", "hairpin_loop", "unpaired_threeprime"),
      count = c(str_count(.x, "f")/total, str_count(.x, "s")/total, str_count(.x, "i")/total, 
                str_count(.x, "m")/total, str_count(.x, "h")/total, str_count(.x, "t")/total)
    )
  }) %>% ldply(rbind)
  
  #make this into a table format that can be easily converted to matrix for heatmap
  table_for_matrix <- tally_per_position %>% 
    pivot_wider(names_from = .id, values_from = count) %>% 
    dplyr::select(-`0`)
  
  return(table_for_matrix)
}

#this function calculates the proportion paired v unpaired at each position for a dotbracket file 
process_dotbracket <- function(rnafold_output){
  
  dotbracket <- process_rnafold_output(rnafold_output, "db")
  
  db_tally_per_position <- 
    dotbracket %>%
    separate(col = 1, into = as.character(seq(from = 0, to = 200)), "") %>% 
    dplyr::select(-`0`)
  
  #dot is considered 0 and bracket is 1 so it is the proportion of paired 
  db_prop <- 
    db_tally_per_position[,-1] %>% 
    dplyr::mutate_all(~ifelse(. == ".", 0, 1)) %>% 
    dplyr::summarise_all(sum) %>% 
    dplyr::mutate_all(~./nrow(db_tally_per_position))
  
  return(db_prop)
  
}


