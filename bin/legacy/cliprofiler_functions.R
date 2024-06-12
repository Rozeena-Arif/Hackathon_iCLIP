## ============================================================================
## The exonProfile function for GRanges objects.
## ----------------------------------------------------------------------------

#' @import dplyr
#' @import methods
#' @import ggplot2
#' @importFrom GenomicRanges GRanges
#' @importFrom Biostrings getSeq
#' @importFrom BSgenome getBSgenome
#' @importFrom GenomicRanges start

## ============================================================================
## Small functions
## ----------------------------------------------------------------------------

## Length filtering
.lenFilteringMax <- function(anno, maxLength)
{
  anno <- anno[width(anno) <= maxLength]
  return(anno)
}
.lenFilteringMin <- function(anno, minLength)
{
  anno <- anno[width(anno) >= minLength]
  return(anno)
}

## Calculation of exon position
.exonPosition <- function(object)
{
  ##-----calculate the position-----##
  object_p <- object[strand(object) == "+"]
  object_p$exon_map <- (start(object_p) -
                          object_p$exon_S)/object_p$exon_length
  object_n <- object[strand(object) == "-"]
  object_n$exon_map <- (object_n$exon_E -
                          start(object_n))/object_n$exon_length
  
  object <- c(object_p, object_n)
  
  ## Give the no mapped peaks a value 3 for their position
  object$exon_map[object$exon_map == -Inf | object$exon_map == Inf] <- 3
  
  return(object)
}

## Exclude the first and last exon for each transcript
.exonExtract <- function(anno)
{
  anno <- as.data.frame(anno) %>% group_by(transcript_id) %>%
    mutate(maxExon = max(exon_number))
  anno <- makeGRangesFromDataFrame(anno, keep.extra.columns = TRUE)
  ## exclude the first and last exon
  anno <- anno[anno$exon_number != 1 & anno$exon_number != anno$maxExon]
  return(anno)
}

## plot
.exonPlot <- function(df, title)
{
  p1 <- ggplot(df, aes(x = exon_map)) +
    geom_density(adjust = 0.2, color = "Orange") +
    theme_bw() + xlab("Exon") +
    scale_x_continuous(breaks = c(0,1),
                       labels = c("3'SS", "5'SS")) +
    ggtitle(title) + ylab("Density of Peaks")
  return(p1)
}

.exonPlotGroup <- function(df, title)
{
  p1 <- ggplot(df, aes(x = exon_map, color = groupIn)) +
    geom_density(adjust = 0.2) +
    theme_bw() + xlab("Exon") +
    scale_x_continuous(breaks = c(0,1),
                       labels = c("3'SS", "5'SS")) +
    ggtitle(title) + ylab("Density of Peaks")
  return(p1)
}
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## The "exonProfile" methods for GRanges objects.
##

#' @title exonProfile for the GRanges objects
#'
#' @description An function to check the position of peaks in the exonic region.
#'
#' @author You Zhou, Kathi Zarnack
#'
#' @param object A GRanges object which should contains all the peaks that you
#'                want to check
#' @param annotation A path way to the annotation file. The format of the
#'                   annotation file should be gff3 and downloaded from
#'                   https://www.gencodegenes.org/
#' @param title The main title for the output meta gene profile plot.
#' @param group The column name which contains the information of grouping
#'     for making the comparison plot. NA means all the peaks belongs to
#'     the same catagory.
#' @param exlevel A parameter for the annotation filtering. exlevel represents
#'     the level that you would like to exclude. NA means no level filtering
#'     for the annotation file. The level from the annotations refers to
#'     how reliable this annotation is. For more information about level
#'     please check
#'     https://www.gencodegenes.org/pages/data_format.html.
#' @param extranscript_support_level A parameter for the annotation filtering.
#'     extranscript_support_level represents the transcript_support_level
#'     that you would like to exclude (e.g. 4 and 5). NA means no
#'     transcript_support_level filtering for the annotation file.
#'     Transcripts are scored according to how well mRNA and EST alignments
#'     match over its full length. Here the number 6 means the
#'     transcript_support_level NA. For more information about level please
#'     check
#'     https://www.gencodegenes.org/pages/data_format.html.
#' @param maxLength A numeric value which indicate the maximum value of exon
#'                  length for the annotation filtering. Or a NA which will
#'                  turn off the max length annotation filtering.
#' @param minLength A numeric value which indicate the minimum value of exon
#'                  length for the annotation filtering. Or a NA which will
#'                  turn off the min length annotation filtering.
#' @param nomap A logical vector (TRUE or FALSE). It indicates whether you
#'              would like to exclude peaks that cannot assign to annotations
#'              in the plot.
#' @details
#' \itemize{
#'     Here is an explanation of output meta data in the \code{list 1}:
#'     \item \code{center}: The center position of each peaks. This center
#'     position is used for calculating the position of peaks within the
#'     genomic regions.
#'     \item \code{exon_S} and \code{exon_E}: The location of 5' and 3'
#'     splice sites (SS) of the exon.
#'     \item \code{exon_length}: The length of the exon that peak assigned.
#'     \item \code{exon_transcript_id}: The transcript ID for the exon
#'     \item \code{exon_map}: The relative position of each peak. This value
#'     close to 0 means this peak located close to the 3' SS. The position
#'     value close to one means the peak close to the 5' SS. Value 3 means this
#'     peaks can not map to any annotation.
#' }
#'
#' @return A list object, the list 1 contains the information of the
#'         assignment of the peaks and their position value within the exon.
#'         The value close to 1 means the peak close to the 5' splice site.
#'         The list 2 includes the plot of exonProfile.
#' @examples
#' ## Load the test data and get the path to the test gff3 file
#' testpath <- system.file("extdata", package = "cliProfiler")
#' test <- readRDS(file.path(testpath, "test.rds"))
#' test_gff3 <- file.path(testpath, "annotation_test.gff3")
#'
#' output <- exonProfile(test, test_gff3)
#' @export
#'

exonProfile <- function(object, annotation, title="Exon Profile", group=NA,
                        exlevel=NA, extranscript_support_level=NA, maxLength=NA, minLength=NA,
                        nomap=FALSE)
{
  if (missing(object))
    stop("The input GRanges object is missing.")
  if (missing(annotation))
    stop("The path to the gff3 annotation file is missing.")
  if (!isS4(object))
    stop("The input object should be a GRanges object.")
  if (!is.logical(nomap))
    stop("The nomap should be a logical vector (TRUE or FALSE).")
  level <- c(1,2,3,NA)
  tsl <- c(1,2,3,4,5,6,NA)
  if (sum(!exlevel %in% level) > 0 |
      sum(!extranscript_support_level %in% tsl) > 0)
    warning("The exlevel should be a vector includes the value of 1, 2,
                3 or NA. extranscript_support_level should be a vector includes
                value 1, 2, 3, 4, 5, 6 or NA.")
  if (!is.numeric(maxLength)&!is.na(maxLength))
    stop("The maxLength should be a numeric value which indicate the
                maximum length of the exon or NA.")
  if (!is.numeric(minLength)&!is.na(minLength))
    stop("The minLength should be a numeric value which indicate the
                minimum length of the exon or NA.")
  if (!is.na(group) &
      (group %in% colnames(elementMetadata(object)) == FALSE))
    stop("When the option group is used, please make sure that
            the group should be a column name of input GRanges object which
            contains the information of which group this peaks belongs to
            and will be used for seperate the for generating the meta gene
            profile curve for different groups")
  
  ## Center peaks and getting exons
  anno <- annotation
  anno_exon <- anno[anno$type == "exon"]
  anno_exon <- as.data.frame(anno_exon) %>%
    group_by(transcript_id) %>%
    mutate(transcript_length = sum(width)) %>%
    makeGRangesFromDataFrame(., keep.extra.columns = TRUE)
  
  ## Give 6 to the NAs in the transcript support level
  anno_exon$transcript_support_level[
    is.na(anno_exon$transcript_support_level)] <- 6
  anno_exon$transcript_support_level[
    anno_exon$transcript_support_level == "NA"] <- 6
  
  object <- .centerPeaks(object)
  
  ## Annotation filtering
  if (sum(is.na(exlevel)) == 0 ) {
    anno_exon <- .annoFilterLevel(anno_exon, exlevel)
  }
  if (sum(is.na(extranscript_support_level)) == 0) {
    anno_exon <- .annoFilterTSL(anno_exon, extranscript_support_level)
  }
  if (!is.na(maxLength)) {
    anno_exon <- .lenFilteringMax(anno_exon, maxLength)
  }
  if (!is.na(minLength)) {
    anno_exon <- .lenFilteringMin(anno_exon, minLength)
  }
  
  anno_exon <- .exonExtract(anno_exon)
  
  anno_exon <- anno_exon[order(anno_exon$level,
                               anno_exon$transcript_support_level,
                               -anno_exon$transcript_length)]
  o <- findOverlaps(object, anno_exon, select = "first")
  
  object$exon_S <-  start(anno_exon)[o]
  object$exon_E <-  end(anno_exon)[o]
  object$exon_length <-  width(anno_exon)[o]
  object$exon_transcript_id <- anno_exon$transcript_id[o]
  
  object$exon_S[is.na(object$exon_S)] <- 0
  object$exon_E[is.na(object$exon_E)] <- 0
  object$exon_length[is.na(object$exon_length)] <- 0
  object$exon_transcript_id[is.na(object$exon_transcript_id)]<-"NO"
  
  object <- .exonPosition(object)
  
  ## shift back to the original peak
  GenomicRanges::start(object) <- object$oriStart
  GenomicRanges::end(object) <- object$oriEnd
  object$oriStart <- NULL
  object$oriEnd <- NULL
  df <- as.data.frame(object)
  
  if (is.na(group)) {
    if (nomap == FALSE) {
      df <- df[df$exon_map != 3,]
      p1 <- .exonPlot(df, title) + coord_cartesian(xlim = c(0,1))
    }
    if (nomap == TRUE){
      p1 <- .exonPlot(df, title) + coord_cartesian(xlim = c(0,3))
    }
  }
  if (!is.na(group)) {
    df$groupIn <- df[,colnames(df) == group]
    if (nomap == FALSE) {
      df <- df[df$exon_map != 3,]
      p1 <- .exonPlotGroup(df, title) + coord_cartesian(xlim = c(0,1))
    }
    if (nomap == TRUE){
      p1 <- .exonPlotGroup(df, title) + coord_cartesian(xlim = c(0,3))
    }
  }
  
  output <- list(Peaks = object, Plot = p1)
  return(output)
}

## ============================================================================
## The metaGeneProfile function for GRanges objects.
## ----------------------------------------------------------------------------

#' @import dplyr
#' @import methods
#' @import ggplot2
#' @importFrom GenomicRanges width
#' @importFrom GenomicRanges seqnames
#' @importFrom GenomicRanges start
#' @importFrom GenomicRanges strand
#' @importFrom GenomicRanges end
#' @importFrom GenomicRanges findOverlaps
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom GenomicRanges GRanges
#' @importFrom rtracklayer import.gff3
#' @importFrom S4Vectors elementMetadata
#' @importFrom utils globalVariables


## To fix the global variable note
utils::globalVariables(c("transcript_id", "geneType", "Fraction", "Number",
                         "groupIn", "exon_map", "Intron_map", "..count..", "exlevel", "maxLength",
                         "extranscript_support_level", "level","extranscript_support_level",
                         "minLength", "window_map","exon_number", "location", "title","."))

## Extract the center position of all input peaks
.centerPeaks <- function(granges)
{
  ## Store the original information for the peaks
  granges$oriStart <- start(granges)
  granges$oriEnd <- end(granges)
  
  ## Extract the center of all peaks
  granges$half_length <- width(granges)/2
  ## For the peaks which width = 1, keep them as what they are
  granges$half_length[granges$half_length == 0.5] <- 0
  granges$half_length <- round(granges$half_length)
  GenomicRanges::start(granges) <- GenomicRanges::start(granges) +
    granges$half_length
  GenomicRanges::end(granges) <- start(granges)
  granges$half_length <- NULL
  
  granges$center <- start(granges)
  return(granges)
}

## Function for filtering annotations
.annoFilterLevel <- function(anno, exlevel = exlevel)
{
  anno <- anno[!anno$level %in% exlevel]
  return(anno)
}

.annoFilterTSL <- function(anno, extsl = extranscript_support_level)
{
  anno <- anno[!anno$transcript_support_level %in% extsl]
  return(anno)
}

## Extract the position of the each annotation fragment (positive strand)
## without the introns.
.annoCalculateP <- function(annoP){
  annoP <- as.data.frame(annoP)
  annoP <- annoP %>% group_by(transcript_id) %>%
    mutate(full_length = sum(width)) %>%
    mutate(Rstart = min(start)) %>% mutate(Rend = max(end))
  annoP <- makeGRangesFromDataFrame(annoP,keep.extra.columns = TRUE)
  annoP <- sort(annoP) %>% as.data.frame()
  annoP <- annoP %>% group_by(transcript_id) %>%
    mutate(Add = c(0, cumsum(width)[-length(width)]))
  annoP <- makeGRangesFromDataFrame(annoP, keep.extra.columns = TRUE)
  return(annoP)
}

## Extract the position of the each annotation fragment (negative strand)
## without the introns.
.annoCalculateN <- function(annoN)
{
  annoN <- as.data.frame(annoN)
  annoN <- annoN %>% group_by(transcript_id) %>%
    mutate(full_length = sum(width)) %>%
    mutate(Rstart = min(start)) %>% mutate(Rend = max(end))
  annoN <- makeGRangesFromDataFrame(annoN,keep.extra.columns = TRUE)
  annoN <- sort(annoN, decreasing=TRUE) %>% as.data.frame()
  annoN <- annoN %>% group_by(transcript_id) %>%
    mutate(Add = c(0, cumsum(width)[-length(width)]))
  annoN <- makeGRangesFromDataFrame(annoN,keep.extra.columns = TRUE)
  return(annoN)
}

## The functions for assigning peaks
.peakAssignment <- function(granges, overlap, fullAnno)
{
  granges$location <- fullAnno$type2[overlap]
  granges$length <- fullAnno$full_length[overlap]
  granges$Rstart <- fullAnno$Rstart[overlap]
  granges$Rend <- fullAnno$Rend[overlap]
  granges$S1 <- start(fullAnno)[overlap]
  granges$E1 <- end(fullAnno)[overlap]
  granges$Add <- fullAnno$Add[overlap]
  granges$Gene_ID <- fullAnno$gene_id[overlap]
  granges$Transcript_ID <- fullAnno$transcript_id[overlap]
  return(granges)
}

## Calculate the position of peaks without intron
.intronOutPosition <- function(granges)
{
  granges$Position <- 5
  granges_n <- granges[strand(granges) == "-"]
  granges_p <- granges[strand(granges) == "+"]
  
  granges_n$Position <- (granges_n$E1-start(granges_n) +
                           granges_n$Add)/granges_n$length
  granges_p$Position <- (start(granges_p) - granges_p$S1 +
                           granges_p$Add) /granges_p$length
  granges <- c(granges_n,granges_p)
  
  granges$Position[granges$Position == -Inf] <- 5
  granges$Position[granges$Position ==  Inf] <- 5
  granges$length <- NULL
  granges$S1 <- NULL
  granges$E1 <- NULL
  granges$Add <- NULL
  granges$Rstart <- NULL
  granges$Rend <- NULL
  
  return(granges)
}

## Making new annotations for including the intron region
.makeTranscriptAnno <- function(anno)
{
  anno <- as.data.frame(anno)
  anno <- anno %>% group_by(transcript_id) %>%
    mutate(Rstart = min(start)) %>% mutate(Rend = max(end))
  anno$order <- seq_len(nrow(anno))
  anno <- anno %>% group_by(transcript_id) %>% filter(order == min(order))
  anno <- makeGRangesFromDataFrame(anno,keep.extra.columns = TRUE)
  GenomicRanges::start(anno) <- anno$Rstart
  GenomicRanges::end(anno) <- anno$Rend
  return(anno)
}

## Calculate position of peaks for the full transcript region (include intron)
.intronInPosition <- function(granges)
{
  granges$Position <- 5
  granges$all_length <- abs(granges$Rstart-granges$Rend)
  granges_n <- granges[strand(granges) == "-"]
  granges_p <- granges[strand(granges) == "+"]
  
  granges_n$Position <- (granges_n$Rend-end(granges_n))/granges_n$all_length
  granges_p$Position <- (start(granges_p) -
                           granges_p$Rstart)/granges_p$all_length
  granges <- c(granges_n,granges_p)
  granges$Position[granges$Position == -Inf] <- 5
  granges$Position[granges$Position ==  Inf] <- 5
  granges$length <- NULL
  granges$all_length <- NULL
  granges$Rstart <- NULL
  granges$Rend <- NULL
  granges$S1 <- NULL
  granges$E1 <- NULL
  
  return(granges)
}

.plotMeta <- function(df, title, adjust)
{
  df$Position[df$location == "CDS"] <-
    df$Position[df$location == "CDS"] + 1
  df$Position[df$location == "UTR3"] <-
    df$Position[df$location == "UTR3"] + 2
  p1 <- ggplot(df, aes(x = Position)) +
    geom_density(fill = "white", alpha= 0.05,
                 adjust = adjust, color = "Orange")  +
    ggtitle(title) +
    scale_x_continuous(breaks = c(0.5,1.5,2.5, 5),
                       labels = c("5'UTR","CDS","3'UTR", "No Map")) +
    geom_vline(xintercept=c(1,2), linetype=2, color="cornflowerblue") +
    theme_bw()
  return(p1)
}

.plotMetaSplit <- function(df, title, adjust)
{
  df$Position[df$location == "CDS"] <-
    df$Position[df$location == "CDS"] + 1
  df$Position[df$location == "UTR3"] <-
    df$Position[df$location == "UTR3"] + 2
  p1 <- ggplot(df, aes(x = Position, color = location)) +
    geom_density(fill = "white", alpha= 0.05,
                 adjust = adjust)  +
    ggtitle(title) +
    scale_x_continuous(breaks = c(0.5,1.5,2.5, 5),
                       labels = c("5'UTR","CDS","3'UTR", "No Map")) +
    geom_vline(xintercept=c(1,2), linetype=2, color="cornflowerblue") +
    theme_bw() +
    geom_density(data=df, aes(x = Position), color="grey",
                 linetype = 2, fill="grey", alpha=0.1)
  return(p1)
}

.plotMetaGroup <- function(df, title, adjust)
{
  df$Position[df$location == "CDS"] <-
    df$Position[df$location == "CDS"] + 1
  df$Position[df$location == "UTR3"] <-
    df$Position[df$location == "UTR3"] + 2
  p1 <- ggplot(df, aes(x=Position, color=groupIn)) +
    geom_density(fill="white",alpha= 0.05,adjust=adjust)  +
    ggtitle(title) +
    scale_x_continuous(breaks=c(0.5,1.5,2.5,5),
                       labels = c("5'UTR","CDS","3'UTR","No Map"))  +
    geom_vline(xintercept = c(1,2), linetype = 2,
               color = "cornflowerblue") +
    theme_bw() + theme(legend.position="bottom") +
    guides(color=guide_legend(title="Group"))
  return(p1)
}

.plotTranscript <- function(df, title, adjust)
{
  df$Position[df$location == "CDS"] <-
    df$Position[df$location == "CDS"] + 1
  df$Position[df$location == "UTR3"] <-
    df$Position[df$location == "UTR3"] + 2
  p1 <- ggplot(df, aes(x=Position)) +
    geom_density(fill="white", alpha=0.05,
                 adjust = adjust, color="Orange")  +
    ggtitle(title) +
    scale_x_continuous(breaks = c(0,1,2,3, 5),
                       labels = c("5' End","Start Codon", "Stop Codon", "3' End", "No Map")) +
    geom_vline(xintercept=c(1,2), linetype=2, color="cornflowerblue") +
    theme_bw()
  return(p1)
}

.plotTranscriptSplit <- function(df, title, adjust)
{
  df$Position[df$location == "CDS"] <-
    df$Position[df$location == "CDS"] + 1
  df$Position[df$location == "UTR3"] <-
    df$Position[df$location == "UTR3"] + 2
  p1 <- ggplot(df, aes(x = Position, color = location)) +
    geom_density(fill = "white", alpha= 0.05,
                 adjust = adjust)  +
    ggtitle(title) +
    scale_x_continuous(breaks = c(0,1,2,3, 5),
                       labels = c("5' End","Start Codon",
                                  "Stop Codon", "3' End", "No Map")) +
    geom_vline(xintercept=c(1,2), linetype=2, color="cornflowerblue") +
    theme_bw() +
    geom_density(data=df, aes(x = Position), color="grey",
                 linetype = 2, fill="grey", alpha=0.1)
  return(p1)
}

.plotTranscriptGroup <- function(df, title, adjust)
{
  df$Position[df$location == "CDS"] <-
    df$Position[df$location == "CDS"] + 1
  df$Position[df$location == "UTR3"] <-
    df$Position[df$location == "UTR3"] + 2
  p1 <- ggplot(df, aes(x = Position, color = groupIn)) +
    geom_density(fill = "white", alpha= 0.05, adjust = adjust) +
    ggtitle(title) +
    scale_x_continuous(breaks = c(0, 1 ,2 ,3 ,5),
                       labels = c("5' End","Start Codon", "Stop Codon", "3' End", "No Map")) +
    geom_vline(xintercept = c(1,2), linetype = 2,
               color = "cornflowerblue") +
    theme_bw() + theme(legend.position="bottom") +
    guides(color=guide_legend(title="Group"))
  return(p1)
}



## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## The "metaGeneProfile" methods for GRanges objects.
##

#' @title metaGeneProfile for the GRanges objects
#'
#' @description An function for calculating the genomic position and generate
#'              the meta gene profile plot of the input peaks.
#'
#' @author You Zhou, Kathi Zarnack
#'
#' @param object A GRanges object which should contains all the peaks that you
#'                want to check
#' @param annotation A path way to the annotation file. The format of the
#'                   annotation file should be gff3 and downloaded from
#'                   https://www.gencodegenes.org/
#' @param include_intron A logical vector TRUE or FALSE that define whether
#'                       the intronic region should be included in the position
#'                       calculation or not.
#' @param title The main title for the output meta gene profile plot.
#' @param group The column name which contains the information of grouping
#'     for making the comparison plot. NA means all the peaks belongs to
#'     the same catagory.
#' @param split A logical vector which indicates whether the plot should show
#'     the density curve for 3'UTR, CDS 5'UTR, respectively.
#' @param exlevel A parameter for the annotation filtering. exlevel represents
#'     the level that you would like to exclude. NA means no level filtering
#'     for the annotation file. The level from the annotations refers to
#'     how reliable this annotation is. For more information about level
#'     please check
#'     https://www.gencodegenes.org/pages/data_format.html.
#' @param extranscript_support_level A parameter for the annotation filtering.
#'     extranscript_support_level represents the transcript_support_level
#'     that you would like to exclude (e.g. 4 and 5). NA means no
#'     transcript_support_level filtering for the annotation file.
#'     Transcripts are scored according to how well mRNA and EST alignments
#'     match over its full length. Here the number 6 means the
#'     transcript_support_level NA. For more information about level please
#'     check
#'     https://www.gencodegenes.org/pages/data_format.html.
#' @param adjust A parameter inherit from ggplot2. A multiplicate bandwidth
#'               adjustment. This makes it possible to adjust the bandwidth
#'               while still using the a bandwidth estimator. For example,
#'               adjust = 1/2 means use half of the default bandwidth.
#' @param nomap A logical vector. It indicates whether you would like to
#'              exclude peaks that cannot assign to annotations in the plot.
#' @details
#' \itemize{
#'     Here is an explanation of output meta data in the \code{list 1}:
#'     \item \code{center}: The center position of each peaks. This center
#'     position is used for calculating the position of peaks within the
#'     genomic regions.
#'     \item \code{location}: Which genomic region this peak belongs to.
#'     \item \code{Gene ID}: Which gene this peak belongs to.
#'     \item \code{Position}: The relative position of each peak. This
#'     value close to 0 means this peak located close to the 5' end of the
#'     genomic feature. The position value close to one means the peak close to
#'     the 3' end of the genomic feature. Value 5 means this peaks can not map
#'     to any annotation.
#' }
#'
#' @return A list object, the list 1 contains the information of the assignment
#'     of the peaks and their position value. The position value between 0 to 1
#'     means it located at the 5' UTR, the value close to the 1 means the
#'     position of this peak close to the 3' end of the 5' UTR. Peaks located
#'     at CDS would have a number between 1 and 2. Postion value between 2 to 3
#'     means this peak assigned to the 3' UTR. For the peaks which can not be
#'     assignment to any annotations, they have the value 5. The list 2
#'     includes the plot of meta gene profile.
#' @examples
#' ## Load the test data and get the path to the test gff3 file
#' testpath <- system.file("extdata", package = "cliProfiler")
#' test <- readRDS(file.path(testpath, "test.rds"))
#' test_gff3 <- file.path(testpath, "annotation_test.gff3")
#'
#' output <- metaGeneProfile(
#'   object = test, annotation = test_gff3,
#'   include_intron = FALSE
#' )
#' @export
#'

metaGeneProfile <- function(object, annotation, include_intron=FALSE,
                            title="Meta Gene Profile", group=NA, split=FALSE,
                            exlevel=NA, extranscript_support_level=NA,
                            adjust=1, nomap=FALSE)
{
  if(missing(object))
    stop("The input GRanges object is missing.")
  if(!isS4(object))
    stop("The input should be a S4 GRanges object.")
  if(missing(annotation))
    stop("The pathway to the annotation file is needed.")
  if(!is.logical(include_intron))
    stop("The include_intron should be a logical vector,
            i.e. TRUE or FALSE.")
  if(!is.logical(nomap))
    stop("The nomap should be a logical vector, i.e. TRUE or FALSE.")
  if(!is.numeric(adjust))
    stop("The adjust should be a numberic, e.g. 1 or 0.5.")
  if(!is.character(title))
    stop("The title should be a character string.")
  if(length(GenomicRanges::seqnames(object)) == 0)
    stop("The input object should not be empty.")
  if (!is.na(group) &
      (group %in% colnames(elementMetadata(object)) == FALSE))
    stop("When the option group is used, please make sure that
            the group should be a column name of input GRanges object which
            contains the information of which group this peaks belongs to
            and will be used for seperate the for generating the meta gene
            profile curve for different groups")
  level <- c(1,2,3,NA)
  tsl <- c(1,2,3,4,5,6,NA)
  if (sum(!exlevel %in% level) > 0 |
      sum(!extranscript_support_level %in% tsl) > 0)
    warning("The exlevel should be a vector includes the value of 1, 2,
                3 or NA. extranscript_support_level should be a vector includes
                value 1, 2, 3, 4, 5, 6 or NA.")
  if (isTRUE(split) & !is.na(group))
    stop("If you would like to use the group option, please set split to
            FALSE. Conversely, when you set split to TRUE please set group
            to NA to make sure the plot is readable.")
  if(isS4(object) & length(GenomicRanges::seqnames(object)) != 0){
    ## import annotation file
    anno <- annotation 
    
    ## Annotation filtering
    anno$transcript_support_level[
      is.na(anno$transcript_support_level)] <- 6
    anno$transcript_support_level[
      anno$transcript_support_level == "NA"] <- 6
    
    if (sum(is.na(exlevel)) == 0) {
      anno <- .annoFilterLevel(anno, exlevel)
    }
    if (sum(is.na(extranscript_support_level)) == 0) {
      anno <- .annoFilterTSL(anno, extranscript_support_level)
    }
    
    object <- .centerPeaks(object)
    
    CDS <- anno[anno$type == "CDS"]
    UTR3 <- anno[anno$type == "three_prime_UTR"]
    UTR5 <- anno[anno$type == "five_prime_UTR"]
    
    if(include_intron == FALSE){
      ## Separate different annotation regions for the
      ## calculation
      CDS_p  <- CDS[strand(CDS) == "+"] %>% sort()
      CDS_n  <- CDS[strand(CDS) == "-"] %>% sort()
      
      UTR3_p <- UTR3[strand(UTR3) == "+"] %>% sort()
      UTR3_n <- UTR3[strand(UTR3) == "-"] %>% sort()
      
      UTR5_p <- UTR5[strand(UTR5) == "+"] %>% sort()
      UTR5_n <- UTR5[strand(UTR5) == "-"] %>% sort()
      
      ## Get the location of annotation fragments for intronic
      ## region exclusive version.
      ## 5'UTR
      UTR5_p <- .annoCalculateP(UTR5_p)
      UTR5_n <- .annoCalculateN(UTR5_n)
      UTR5 <- sort(c(UTR5_n,UTR5_p))
      UTR5$type2 <- "UTR5"
      
      ## CDS
      CDS_p <- .annoCalculateP(CDS_p)
      CDS_n <- .annoCalculateN(CDS_n)
      CDS <- sort(c(CDS_n,CDS_p))
      CDS$type2 <- "CDS"
      
      ## 3' UTR
      UTR3_p <- .annoCalculateP(UTR3_p)
      UTR3_n <- .annoCalculateN(UTR3_n)
      UTR3 <- sort(c(UTR3_n,UTR3_p))
      UTR3$type2 <- "UTR3"
      ## Assign the transcript length for the following peak
      ## assignment
      anno_exon <- anno[anno$type == "exon"]
      anno_exon <- as.data.frame(anno_exon) %>%
        group_by(transcript_id) %>%
        mutate(trans_len = sum(width)) %>%
        makeGRangesFromDataFrame(., keep.extra.columns = TRUE)
      
      UTR5$transcript_length <- anno_exon$trans_len[
        match(UTR5$transcript_id, anno_exon$transcript_id)]
      UTR3$transcript_length <- anno_exon$trans_len[
        match(UTR3$transcript_id, anno_exon$transcript_id)]
      CDS$transcript_length <- anno_exon$trans_len[
        match(CDS$transcript_id, anno_exon$transcript_id)]
      
      ## Rank the transcript fragments base on the level,
      ## transcript support level and transcript length ##
      fullAnno <- c(UTR5,UTR3,CDS)
      
      fullAnno <-
        fullAnno[order(fullAnno$level,
                       fullAnno$transcript_support_level,
                       -fullAnno$transcript_length)]
      
      ## Assign peaks to annotation according to
      ## level->transcript->transcript length
      ## Use select = "first" to assign peaks to the first
      ## ranked annotation
      o <- findOverlaps(object, fullAnno, select = "first")
      object <- .peakAssignment(object, o, fullAnno)
      
      ## Remove all the NA from the object
      object$location[is.na(object$location)] <- "NO"
      object$length[is.na(object$length)] <- 0
      object$Rstart[is.na(object$Rstart)] <- 0
      object$Rend[is.na(object$Rend)] <- 0
      object$S1[is.na(object$S1)] <- 0
      object$E1[is.na(object$E1)] <- 0
      object$Add[is.na(object$Add)] <- 0
      object$Gene_ID[is.na(object$Gene_ID)] <- "Nan"
      
      GenomicRanges::start(object) <- object$oriStart
      GenomicRanges::end(object) <- object$oriEnd
      object$oriStart <- NULL
      object$oriEnd <- NULL
      
      ## Get result for the intron excluded version
      object <- .intronOutPosition(object)
      df <- as.data.frame(object)
      if (is.na(group)) {
        if (split == TRUE) {
          if(nomap == FALSE){
            df <- df[df$location != "NO",]
            df$location <- factor(df$location,
                                  levels = c("UTR5", "CDS", "UTR3"))
          }
          ourpic <- .plotMetaSplit(df, title=title, adjust=adjust)
        }
        if (split == FALSE) {
          ourpic <- .plotMeta(df, title=title, adjust=adjust)
        }
      }
      if (!is.na(group)) {
        df$groupIn <- df[,colnames(df) == group]
        if (nomap == FALSE) {
          ourpic <- .plotMetaGroup(df, title=title,
                                   adjust=adjust)
        }
        if (nomap == TRUE) {
          ourpic <- .plotMetaGroup(df, title=title, adjust=adjust)
        }
      }
      if (nomap == FALSE) {
        ourpic <- ourpic + coord_cartesian(xlim = c(0,3))
      }
      if (nomap == TRUE) {
        ourpic <- ourpic + coord_cartesian(xlim = c(0,5))
      }
      ## get output
      output <- list(Peaks = object, Plot = ourpic)
    }
    if (include_intron == TRUE) {
      
      ## Making new annotation includes intronic region for the
      ## peaks assignment
      UTR5 <- .makeTranscriptAnno(UTR5)
      UTR5$type2 <- "UTR5"
      UTR3 <- .makeTranscriptAnno(UTR3)
      UTR3$type2 <- "UTR3"
      CDS  <- .makeTranscriptAnno(CDS)
      CDS$type2 <- "CDS"
      
      anno.transcript <- anno[anno$type == "transcript"]
      UTR5$transcript_length <-
        width(anno.transcript[match(UTR5$transcript_id,
                                    anno.transcript$transcript_id)])
      UTR3$transcript_length <-
        width(anno.transcript[match(UTR3$transcript_id,
                                    anno.transcript$transcript_id)])
      CDS$transcript_length <-
        width(anno.transcript[match(CDS$transcript_id,
                                    anno.transcript$transcript_id)])
      
      ## Making all the annotation together and rank them
      ## according to level->transcript->transcript length
      fullAnno <- c(UTR5,UTR3,CDS)
      
      fullAnno <-
        fullAnno[order(fullAnno$level,
                       fullAnno$transcript_support_level,
                       -fullAnno$transcript_length)]
      
      ## Assign peaks
      o <- findOverlaps(object, fullAnno, select = "first")
      object <- .peakAssignment(object, o, fullAnno)
      
      ## Remove all the NA from the object
      object$location[is.na(object$location)] <- "NO"
      object$Rstart[is.na(object$Rstart)] <- 0
      object$Rend[is.na(object$Rend)] <- 0
      object$S1[is.na(object$S1)] <- 0
      object$E1[is.na(object$E1)] <- 0
      object$Gene_ID[is.na(object$Gene_ID)] <- "Nan"
      
      object <- .intronInPosition(object)
      
      GenomicRanges::start(object) <- object$oriStart
      GenomicRanges::end(object) <- object$oriEnd
      object$oriStart <- NULL
      object$oriEnd <- NULL
      
      df <- as.data.frame(object)
      if (is.na(group)) {
        if (split == TRUE) {
          if(nomap == FALSE){
            df <- df[df$location != "NO",]
            df$location <- factor(df$location,
                                  levels = c("UTR5", "CDS", "UTR3"))
          }
          ourpic <- .plotTranscriptSplit(df, title=title,
                                         adjust=adjust)
        }
        if (split == FALSE) {
          ourpic <- .plotTranscript(df, title=title,
                                    adjust=adjust)
        }
      }
      if (!is.na(group)) {
        df$groupIn <- df[,colnames(df) == group]
        ourpic <- .plotTranscriptGroup(df, title=title,
                                       adjust=adjust)
      }
      if (nomap == FALSE) {
        ourpic <- ourpic + coord_cartesian(xlim = c(0,3))
      }
      if (nomap == TRUE) {
        ourpic <- ourpic + coord_cartesian(xlim = c(0,5))
      }
      output <- list(Peaks = object, Plot = ourpic)
    }
    return(output)
  }
}

## ============================================================================
## The spliceSiteProfile function for GRanges objects.
## ----------------------------------------------------------------------------

#' @import dplyr
#' @import methods
#' @import ggplot2
#' @importFrom GenomicRanges width
#' @importFrom GenomicRanges seqnames
#' @importFrom GenomicRanges start
#' @importFrom GenomicRanges strand
#' @importFrom GenomicRanges end
#' @importFrom GenomicRanges findOverlaps
#' @importFrom GenomicRanges makeGRangesFromDataFrame
#' @importFrom GenomicRanges GRanges
#' @importFrom rtracklayer import.gff3
#' @importFrom S4Vectors elementMetadata

## ============================================================================
## Small functions
## ----------------------------------------------------------------------------

## build the GRange objects of splice site from annotation file
.buildGRend <- function(anno)
{
  anno <- data.frame(seqnames=GenomicRanges::seqnames(anno),
                     start=GenomicRanges::end(anno),
                     end=GenomicRanges::end(anno),
                     strand=GenomicRanges::strand(anno),
                     transcript_support_level=anno$transcript_support_level,
                     level=anno$level,
                     transcript_length=anno$transcript_length,
                     ss=GenomicRanges::end(anno))
  anno <- makeGRangesFromDataFrame(anno,keep.extra.columns=TRUE)
  return(anno)
}

.buildGRstart <- function(anno)
{
  anno <- data.frame(seqnames=GenomicRanges::seqnames(anno),
                     start=GenomicRanges::start(anno),
                     end=GenomicRanges::start(anno),
                     strand=GenomicRanges::strand(anno),
                     transcript_support_level=anno$transcript_support_level,
                     level=anno$level,
                     transcript_length=anno$transcript_length,
                     ss=GenomicRanges::start(anno))
  anno <- makeGRangesFromDataFrame(anno, keep.extra.columns=TRUE)
  return(anno)
}
## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## The "spliceSiteProfile" methods for GRanges objects.
##

#' @title spliceSiteProfile for the GRanges objects
#'
#' @description An function to check the enrichment of peaks around the splice
#'              sites in a absolute distance.
#'
#' @author You Zhou, Kathi Zarnack
#'
#' @param object A GRanges object which should contains all the peaks that you
#'                want to check
#' @param annotation A path way to the annotation file. The format of the
#'                   annotation file should be gff3 and downloaded from
#'                   https://www.gencodegenes.org/
#' @param title The main title for the output meta gene profile plot.
#' @param exlevel A parameter for the annotation filtering. exlevel represents
#'     the level that you would like to exclude. NA means no level filtering
#'     for the annotation file. The level from the annotations refers to
#'     how reliable this annotation is. For more information about level
#'     please check
#'     https://www.gencodegenes.org/pages/data_format.html.
#' @param extranscript_support_level A parameter for the annotation filtering.
#'     extranscript_support_level represents the transcript_support_level
#'     that you would like to exclude (e.g. 4 and 5). NA means no
#'     transcript_support_level filtering for the annotation file.
#'     Transcripts are scored according to how well mRNA and EST alignments
#'     match over its full length. Here the number 6 means the
#'     transcript_support_level NA. For more information about level please
#'     check
#'     https://www.gencodegenes.org/pages/data_format.html.
#' @param exon_length_filtering The exon_length_filtering should be a logical
#'     value which indicated whether user would like to exclude the exons
#'     that have a length less than flanking value. Set this parameter to TRUE
#'     to turn on this filtering step.
#' @param intron_length_filtering The intron_length_filtering should be a
#'     logical value which indicated whether user would like to exclude the
#'     introns that have a length less than flanking value. Set this parameter
#'     to TRUE to turn on this filtering step.
#' @param flanking The size of the flanking windows that you would like to
#'                 check. Flanking=5 will give you the result of the 10+1nt
#'                 windows around the center of peaks.
#' @param bin A number that indicates how many bins would you like to use in
#'            the histogram.
#'
#' @return A list object, the list 1 contains the information of the
#'         position of peaks around 5' or 3' splice sites. The list 2
#'         includes the plot of spliceSiteProfile
#' @examples
#' ## Load the test data and get the path to the test gff3 file
#' testpath <- system.file("extdata", package = "cliProfiler")
#' test <- readRDS(file.path(testpath, "test.rds"))
#' test_gff3 <- file.path(testpath, "annotation_test.gff3")
#'
#' output <- spliceSiteProfile(test, test_gff3,
#'   flanking = 200, bin = 40
#' )
#' @export

spliceSiteProfile <- function(object, annotation, title="Splice Site Profile",
                              exlevel=NA, extranscript_support_level=NA, exon_length_filtering=TRUE,
                              intron_length_filtering=TRUE, flanking=150, bin=30)
{
  if (missing(object))
    stop("The input GRanges object is missing.")
  if (missing(annotation))
    stop("The path to the gff3 annotation file is missing.")
  if (!isS4(object))
    stop("The input object should be a GRanges object.")
  if (!is.numeric(flanking))
    stop("The flanking should be a numeric.")
  if (!is.numeric(bin))
    stop("The parameter bin should be a numeric.")
  level <- c(1,2,3,NA)
  tsl <- c(1,2,3,4,5,6,NA)
  if (sum(!exlevel %in% level) > 0 |
      sum(!extranscript_support_level %in% tsl) > 0)
    warning("The exlevel should be a vector includes the value of 1, 2,
                3 or NA. extranscript_support_level should be a vector includes
                value 1, 2, 3, 4, 5, 6 or NA.")
  if (!is.logical(exon_length_filtering))
    stop("The exon_length_filtering should be a logical value.")
  if (!is.logical(intron_length_filtering))
    stop("The intron_length_filtering should be a logical value.")
  
  ## Center peaks and getting exons
  anno <- annotation
  trans <- anno[anno$type=="transcript"]
  introns <- .getIntrons(anno)
  
  ## Give 6 to the NAs in the transcript support level
  introns$transcript_support_level[
    is.na(introns$transcript_support_level)] <- 6
  introns$transcript_support_level[
    introns$transcript_support_level == "NA"] <- 6
  introns$transcript_length <-
    GenomicRanges::width(trans)[match(introns$transcript_id,
                                      trans$transcript_id)]
  ## center peaks
  object <- .centerPeaks(object)
  ## annotation filtering
  if (sum(is.na(exlevel)) == 0) {
    introns <- .annoFilterLevel(introns, exlevel)
  }
  if (sum(is.na(extranscript_support_level)) == 0) {
    introns <- .annoFilterTSL(introns, extranscript_support_level)
  }
  if (isTRUE(exon_length_filtering)) {
    introns <- introns[introns$exon1_l >= flanking]
    introns <- introns[introns$exon2_l >= flanking]
  }
  if (isTRUE(intron_length_filtering)) {
    introns <- introns[width(introns) >= flanking]
  }
  
  ## split positive and negative strand
  introns_p <- introns[strand(introns) == "+"]
  introns_n <- introns[strand(introns) == "-"]
  
  ## Check the 5' splice site (5'SS)
  ## First generate the coordinates for the 5'SS
  start_p <- .buildGRstart(introns_p)
  start_n <- .buildGRend(introns_n)
  ## Filter out the duplicates
  
  start <- c(start_p, start_n)
  start <- unique(start)
  
  start <- start[order(start$level,
                       start$transcript_support_level,
                       -start$transcript_length)]
  
  ## Extend the 5'SS for checking
  start <- start+flanking
  o <- findOverlaps(object, start, select = "first")
  
  object$start_2 <- start$ss[o]
  
  ## Check 3' splice sites (3'SS)
  end_p <- .buildGRend(introns_p)
  end_n <- .buildGRstart(introns_n)
  end <- c(end_p, end_n)
  end <- unique(end)
  
  end <- end[order(end$level,
                   end$transcript_support_level,
                   -end$transcript_length)]
  
  end <- end+flanking
  oo <- findOverlaps(object, end, select = "first")
  
  object$end_2 <- end$ss[oo]
  
  object_p <- object[strand(object)=="+"]
  object_n <- object[strand(object)=="-"]
  
  ## calculation
  object_p$Position5SS <- start(object_p)-object_p$start_2
  object_n$Position5SS <- object_n$start_2-start(object_n)
  
  object_p$Position3SS <- start(object_p)-object_p$end_2
  object_n$Position3SS <- object_n$end_2-start(object_n)
  
  object <- c(object_n, object_p)
  
  ## Make dataframe without NAs for the ggplot
  df_5SS <- data.frame(Position=
                         object$Position5SS[!is.na(object$Position5SS)],
                       Type="5' Splice Site")
  df_3SS <- data.frame(Position=
                         object$Position3SS[!is.na(object$Position3SS)],
                       Type="3' Splice Site")
  
  df <- rbind(df_5SS, df_3SS)
  
  df$Type <- factor(df$Type, levels = c("5' Splice Site",
                                        "3' Splice Site"))
  p <- ggplot(df, aes(x = Position)) +
    geom_histogram(bins = bin, fill = "orange", color = "grey") +
    theme_bw() + facet_wrap(~Type) + ggtitle(title) +
    geom_vline(xintercept = 0, linetype = 2,
               color = "cornflowerblue")
  
  object$start_2 <- NULL
  object$end_2 <- NULL
  GenomicRanges::start(object) <- object$oriStart
  GenomicRanges::end(object) <- object$oriEnd
  object$oriStart <- NULL
  object$oriEnd <- NULL
  
  out <- list(Peaks=object, Plot=p)
  return(out)
}



## ============================================================================
## The intronProfile function for GRanges objects.
## ----------------------------------------------------------------------------
#' @import dplyr
#' @import methods
#' @import ggplot2
#' @importFrom S4Vectors split
#' @importFrom GenomicRanges psetdiff
#' @importFrom GenomicRanges width
#' @importFrom GenomicRanges start
#' @importFrom GenomicRanges strand
#' @importFrom GenomicRanges end

## ============================================================================
## Small functions
## ----------------------------------------------------------------------------

.getIntrons <- function(anno)
{
  trans <- anno[anno$type=="transcript"]
  exon.gr <- anno[anno$type=="exon"]
  exon.gr$exonID <- paste0(exon.gr$transcript_id, ",", exon.gr$exon_number)
  exonsByTranscript <- split(anno[anno$type=="exon"],
                             anno[anno$type=="exon"]$transcript_id)
  trans <- trans[match(names(exonsByTranscript),
                       trans$transcript_id)]
  
  introns <- psetdiff(trans, exonsByTranscript)
  introns <- as(introns, "GRangesList") %>% unlist
  
  introns$transcript_id <- names(introns)
  names(introns) <- NULL
  
  ## make intron number
  introns.n <- introns[strand(introns) == "-"] %>%
    sort(., decreasing = TRUE) %>% as.data.frame() %>%
    group_by(transcript_id) %>%
    mutate(intron_number = seq_len(n())) %>%
    makeGRangesFromDataFrame(., keep.extra.columns = TRUE)
  introns.p <- introns[strand(introns) == "+"] %>%
    sort(.) %>% as.data.frame() %>% group_by(transcript_id) %>%
    mutate(intron_number = seq_len(n())) %>%
    makeGRangesFromDataFrame(., keep.extra.columns = TRUE)
  introns <- c(introns.n, introns.p)
  
  introns$intronID1 <- paste0(introns$transcript_id,
                              ",",introns$intron_number)
  introns$intronID2 <- paste0(introns$transcript_id,
                              ",",introns$intron_number+1)
  ## Assign the length of flanking exons this value will be used in ssprofile
  introns$exon1_l <- width(exon.gr)[match(introns$intronID1,
                                          exon.gr$exonID)]
  introns$exon2_l <- width(exon.gr)[match(introns$intronID2,
                                          exon.gr$exonID)]
  
  introns$transcript_length <-
    width(trans)[match(introns$transcript_id,
                       trans$transcript_id)]
  introns$level <-
    trans$level[match(introns$transcript_id,
                      trans$transcript_id)]
  
  ## tsl means transcript_support_level
  introns$tsl <-
    trans$transcript_support_level[match(introns$transcript_id,
                                         trans$transcript_id)]
  
  ## Remove NA avoid warning message
  introns$level[is.na(introns$level)] <- 6
  introns$tsl[is.na(introns$tsl)] <- 6
  introns$tsl[introns$tsl == "NA"] <-6
  
  ## Change the type of data for the ranking step
  introns$level <- as.integer(introns$level)
  introns$transcript_support_level <- as.integer(introns$tsl)
  introns$tsl <- NULL
  
  return(introns)
}

.intronPosition <- function(object)
{
  ##-----calculate the position-----##
  object_p <- object[strand(object) == "+"]
  object_p$Intron_map <- (start(object_p) -
                            object_p$Intron_S)/object_p$Intron_length
  object_n <- object[strand(object) == "-"]
  object_n$Intron_map <- (object_n$Intron_E -
                            start(object_n))/object_n$Intron_length
  
  object <- c(object_p, object_n)
  
  ## Give the no mapped peaks a value 3 for their position
  object$Intron_map[object$Intron_map == -Inf|object$Intron_map == Inf] <- 3
  
  return(object)
}

## Plot
.intronPlot <- function(df, title)
{
  p1 <- ggplot(df, aes(x = Intron_map)) +
    geom_density(adjust = 0.2, color = "Orange") +
    theme_bw() + xlab("Intron") +
    scale_x_continuous(breaks = c(0,1),
                       labels = c("5'SS", "3'SS")) +
    ggtitle(title) + ylab("Density of Peaks")
  return(p1)
}

.intronPlotGroup <- function(df, title)
{
  p1 <- ggplot(df, aes(x = Intron_map, color = groupIn)) +
    geom_density(adjust = 0.2) +
    theme_bw() + xlab("Intron") +
    scale_x_continuous(breaks = c(0,1),
                       labels = c("5'SS", "3'SS")) +
    ggtitle(title) + ylab("Density of Peaks")
  return(p1)
}

## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## The "intronProfile" methods for GRanges objects.
##

#' @title intronProfile for the GRanges objects
#'
#' @description An function to check the position of peaks in the intronic
#'              region.
#'
#' @author You Zhou, Kathi Zarnack
#'
#' @param object A GRanges object which should contains all the peaks that you
#'                want to check
#' @param annotation A path way to the annotation file. The format of the
#'                   annotation file should be gff3 and downloaded from
#'                   https://www.gencodegenes.org/
#' @param title The main title for the output meta gene profile plot.
#' @param group The column name which contains the information of grouping
#'     for making the comparison plot. NA means all the peaks belongs to
#'     the same catagory.
#' @param exlevel A parameter for the annotation filtering. exlevel represents
#'     the level that you would like to exclude. NA means no level filtering
#'     for the annotation file. The level from the annotations refers to
#'     how reliable this annotation is. For more information about level
#'     please check
#'     https://www.gencodegenes.org/pages/data_format.html.
#' @param extranscript_support_level A parameter for the annotation filtering.
#'     extranscript_support_level represents the transcript_support_level
#'     that you would like to exclude (e.g. 4 and 5). NA means no
#'     transcript_support_level filtering for the annotation file.
#'     Transcripts are scored according to how well mRNA and EST alignments
#'     match over its full length. Here the number 6 means the
#'     transcript_support_level NA. For more information about level please
#'     check
#'     https://www.gencodegenes.org/pages/data_format.html.
#' @param maxLength A numeric value which indicate the maximum value of exon
#'                  length for the annotation filtering. Or a NA which will
#'                  turn off the max length annotation filtering.
#' @param minLength A numeric value which indicate the minimum value of exon
#'                  length for the annotation filtering. Or a NA which will
#'                  turn off the min length annotation filtering.
#' @param nomap A logical vector (TRUE or FALSE). It indicates whether you
#'              would like to exclude peaks that cannot assign to annotations
#'              in the plot.
#' @details
#' \itemize{
#'     Here is an explanation of output meta data in the \code{list 1}:
#'     \item \code{center}: The center position of each peaks. This center
#'     position is used for calculating the position of peaks within the
#'     genomic regions.
#'     \item \code{Intron_S} and \code{Intron_E}: The location of 5' and 3'
#'     splice sites (SS) of the intron.
#'     \item \code{Intron_length}: The length of the intron that peak assigned.
#'     \item \code{Intron_transcript_id}: The transcript ID for the intron.
#'     \item \code{Intron_map}: The relative position of each peak. This value
#'     close to 0 means this peak located close to the 5' SS. The position
#'     value close to one means the peak close to the 3' SS. Value 3 means
#'     this peaks can not map to any annotation.
#' }
#'
#' @return A list object, the list 1 contains the information of the
#'     assignment of the peaks and their position value within the intron.
#'     The value close to 1 means the peak close to the 3' splice site.
#'     The list 2 includes the plot of intronProfile.
#' @examples
#' ## Load the test data and get the path to the test gff3 file
#' testpath <- system.file("extdata", package = "cliProfiler")
#' test <- readRDS(file.path(testpath, "test.rds"))
#' test_gff3 <- file.path(testpath, "annotation_test.gff3")
#'
#' output <- intronProfile(test, test_gff3)
#' @export
#'

intronProfile <- function(object, annotation, title="Intron Profile", group=NA,
                          exlevel=NA, extranscript_support_level=NA, maxLength=NA, minLength=NA,
                          nomap=FALSE)
{
  if (missing(object))
    stop("The input GRanges object is missing.")
  if (missing(annotation))
    stop("The path to the gff3 annotation file is missing.")
  if (!isS4(object))
    stop("The input object should be a GRanges object.")
  if (!is.logical(nomap))
    stop("The nomap should be a logical vector (TRUE or FALSE).")
  level <- c(1,2,3,NA)
  tsl <- c(1,2,3,4,5,6,NA)
  if (sum(!exlevel %in% level) > 0 |
      sum(!extranscript_support_level %in% tsl) > 0)
    warning("The exlevel should be a vector includes the value of 1, 2,
                3 or NA. extranscript_support_level should be a vector includes
                value 1, 2, 3, 4, 5, 6 or NA.")
  if (!is.numeric(maxLength)&!is.na(maxLength))
    stop("The maxLength should be a numeric value which indicate the
                maximum length of the exon or NA.")
  if (!is.numeric(minLength)&!is.na(minLength))
    stop("The minLength should be a numeric value which indicate the
                minimum length of the exon or NA.")
  if (!is.na(group) &
      (group %in% colnames(elementMetadata(object)) == FALSE))
    stop("When the option group is used, please make sure that
            the group should be a column name of input GRanges object which
            contains the information of which group this peaks belongs to
            and will be used for seperate the for generating the meta gene
            profile curve for different groups")
  
  ## Center peaks and getting introns
  anno <- rtracklayer::import.gff3(con = annotation)
  object <- .centerPeaks(object)
  introns <- .getIntrons(anno)
  introns$intronID <- seq_len(length(introns))
  
  ## Annotation filtering
  if (sum(is.na(exlevel)) == 0) {
    introns <- .annoFilterLevel(introns, exlevel)
  }
  if (sum(is.na(extranscript_support_level)) == 0) {
    introns <- .annoFilterTSL(introns, extranscript_support_level)
  }
  if (!is.na(maxLength)) {
    introns <- .lenFilteringMax(introns, maxLength)
  }
  if (!is.na(minLength)) {
    introns <- .lenFilteringMin(introns, minLength)
  }
  
  ## Ranking intron for the peaks assignment
  introns <- introns[order(introns$level,
                           introns$transcript_support_level, -introns$transcript_length)]
  
  o2 <- findOverlaps(object, introns, select = "first")
  
  ## Assign the annotation details to the input GRanges
  object$Intron_S <-  start(introns)[o2]
  object$Intron_E <-  end(introns)[o2]
  object$Intron_length <-  width(introns)[o2]
  object$intronID <- introns$intronID[o2]
  object$Intron_transcript_id <- introns$transcript_id[o2]
  
  object$Intron_S[is.na(object$Intron_S)] <- 0
  object$Intron_E[is.na(object$Intron_E)] <- 0
  object$Intron_length[is.na(object$Intron_length)] <- 0
  object$intronID[is.na(object$intronID)] <- "NO"
  object$Intron_transcript_id[is.na(object$Intron_transcript_id)]<-
    "NO"
  
  object <- .intronPosition(object)
  object$intronID <- NULL
  
  GenomicRanges::start(object) <- object$oriStart
  GenomicRanges::end(object) <- object$oriEnd
  object$oriStart <- NULL
  object$oriEnd <- NULL
  
  df <- as.data.frame(object)
  if (is.na(group)) {
    if (nomap == FALSE) {
      df <- df[df$Intron_map != 3,]
      p1 <- .intronPlot(df, title) +
        coord_cartesian(xlim = c(0,1))
    }
    if (nomap == TRUE){
      p1 <- .intronPlot(df, title) +
        coord_cartesian(xlim = c(0,3))
    }
  }
  if (!is.na(group)) {
    df$groupIn <- df[,colnames(df) == group]
    if (nomap == FALSE) {
      df <- df[df$Intron_map != 3,]
      p1 <- .intronPlotGroup(df, title) +
        coord_cartesian(xlim = c(0,1))
    }
    if (nomap == TRUE){
      p1 <- .intronPlotGroup(df, title) +
        coord_cartesian(xlim = c(0,3))
    }
  }
  
  output <- list(Peaks = object, Plot = p1)
  return(output)
}

## ============================================================================
## The geneTypeProfile function for GRanges objects.
## ----------------------------------------------------------------------------

#' @import dplyr
#' @import methods
#' @import ggplot2
#' @importFrom GenomicRanges GRanges


## ============================================================================
## Small functions
## ----------------------------------------------------------------------------

.geneTypePlot <- function(df, title)
{
  p1 <- ggplot(df, aes(x=geneType)) +
    geom_bar(fill="orange", color="cornflowerblue") +
    geom_label(stat='count', aes(label=..count..), vjust=0.8) +
    theme_bw() + ggtitle(title) +
    theme(axis.text.x = element_text(angle = 45, hjust=1))
  return(p1)
}

## - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
## The "geneTypeProfile" methods for GRanges objects.
##

#' @title geneTypeProfile for the GRanges objects
#'
#' @description An function to check the gene type belonging for the peaks.
#'
#' @author You Zhou, Kathi Zarnack
#'
#' @param object A GRanges object which should contains all the peaks that you
#'                want to check
#' @param annotation A path way to the annotation file. The format of the
#'                   annotation file should be gff3 and downloaded from
#'                   https://www.gencodegenes.org/
#' @param title The main title for the output meta gene profile plot.
#' @param exlevel A parameter for the annotation filtering. exlevel represents
#'     the level that you would like to exclude. NA means no level filtering
#'     for the annotation file. The level from the annotations refers to
#'     how reliable this annotation is. For more information about level
#'     please check
#'     https://www.gencodegenes.org/pages/data_format.html.
#' @param extranscript_support_level A parameter for the annotation filtering.
#'     extranscript_support_level represents the transcript_support_level
#'     that you would like to exclude (e.g. 4 and 5). NA means no
#'     transcript_support_level filtering for the annotation file.
#'     Transcripts are scored according to how well mRNA and EST alignments
#'     match over its full length. Here the number 6 means the
#'     transcript_support_level NA. For more information about level please
#'     check
#'     https://www.gencodegenes.org/pages/data_format.html.
#' @details
#' \itemize{
#'     Here is an explanation of output meta data in the \code{list 1}:
#'     \item \code{center}: The center position of each peaks. This center
#'     position is used for calculating the position of peaks within the
#'     genomic regions.
#'     \item \code{geneType}: The gene type of the gene that input peak belongs
#'     to.
#'     \item \code{Gene_ID}: The gene ID of the gene that input peak
#'     belongs to.
#' }
#'
#' @return A list object, the list 1 contains the information of the
#'         assignment of the peaks and the gene type of their located genes.
#'         The list 2 includes the plot of geneTypeProfile
#' @examples
#' ## Load the test data and get the path to the test gff3 file
#' testpath <- system.file("extdata", package = "cliProfiler")
#' test <- readRDS(file.path(testpath, "test.rds"))
#' test_gff3 <- file.path(testpath, "annotation_test.gff3")
#'
#' output <- geneTypeProfile(test, test_gff3)
#' @export

geneTypeProfile <- function(object, annotation, title="Gene Type Profile",
                            exlevel=NA, extranscript_support_level=NA)
{
  if (missing(object))
    stop("The input GRanges object is missing.")
  if (missing(annotation))
    stop("The path to the gff3 annotation file is missing.")
  if (!isS4(object))
    stop("The input object should be a GRanges object.")
  
  ## Center peaks and getting transcripts
  object <- .centerPeaks(object)
  
  anno <- rtracklayer::import.gff3(con = annotation)
  trans <- anno[anno$type=="transcript"]
  trans$transcript_support_level[
    is.na(trans$transcript_support_level)] <- 6
  trans$transcript_support_level[
    trans$transcript_support_level == "NA"] <- 6
  anno_exon <- anno[anno$type == "exon"]
  anno_exon <- as.data.frame(anno_exon) %>%
    group_by(transcript_id) %>%
    mutate(trans_len = sum(width))  %>%
    filter(exon_number == "1") %>%
    makeGRangesFromDataFrame(., keep.extra.columns = TRUE)
  
  trans$transcript_length <- anno_exon$trans_len[
    match(trans$transcript_id, anno_exon$transcript_id)]
  
  ## annotation filtering
  if (sum(is.na(exlevel)) == 0) {
    trans <- .annoFilterLevel(trans, exlevel)
  }
  if (sum(is.na(extranscript_support_level)) == 0) {
    trans <- .annoFilterTSL(trans, extranscript_support_level)
  }
  
  ## Peaks assignment
  trans <- trans[order(trans$level,
                       trans$transcript_support_level, -trans$transcript_length)]
  
  o <- findOverlaps(object, trans, select = "first")
  object$geneType <- trans$gene_type[o]
  object$Gene_ID <- trans$gene_id[o]
  df <- as.data.frame(object)
  
  p1 <- .geneTypePlot(df, title)
  
  GenomicRanges::start(object) <- object$oriStart
  GenomicRanges::end(object) <- object$oriEnd
  object$oriStart <- NULL
  object$oriEnd <- NULL
  
  out <- list(Peaks = object, Plot = p1)
  return(out)
}





