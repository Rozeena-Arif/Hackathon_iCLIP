
.libPaths("/home4/2711498i/R_packages")

library(tidyverse)
library(data.table)
library(IHW)
library(DEWSeq)
library(ggrepel)
library(BiocParallel)
library(scales)

dir.create("./DEW-seq/Regions_w50s20")

## Read data
print("read data")
countData <- fread("DEW-seq/Matrix_w50s20/Full_PURAB_iCLIP.matrix.txt.gz", sep = "\t")

print("read annotation")
annotationData <- fread("GENOMEDIR/sliding_window/HS.GRCh38.SINV_attribute-fix.flattened.w50s20.maptoid.txt.gz", sep = "\t")

extract_regions <- function(resultWindows, padjThresh, log2FoldChangeThresh, trackName, fileName) {
  resultWindows <- as.data.frame(resultWindows)
  resultRegions <- extractRegions(windowRes  = resultWindows,
                                  padjCol    = "p_adj_IHW",
                                  padjThresh = padjThresh, 
                                  log2FoldChangeThresh = log2FoldChangeThresh)
  toBED(windowRes = resultWindows,
        regionRes = resultRegions,
        padjCol    = "p_adj_IHW",
        padjThresh = padjThresh,
        log2FoldChangeThresh = log2FoldChangeThresh,
        trackName = trackName,
        fileName  = fileName
  )
  return(resultRegions)
}

# * * * * * PURB
print("PURA")
##
print("prep data")
PURA <- countData %>%
  dplyr::select(unique_id, contains("PURA"), -contains("GFP"), -contains("PURB")) %>%
  #dplyr::select(unique_id, contains("IP"), contains("SMI")) %>%
  column_to_rownames(var = "unique_id")

PURA_colData <- tibble(name = colnames(PURA)) %>%
  mutate(type = case_when(
    #str_detect(name, "IP") ~ "IP",
    str_detect(name, "SMI") ~ "SMI",
    TRUE ~ "IP"
  )) %>%
  mutate(type = fct_relevel(type, c("IP", "SMI"))) %>%
  column_to_rownames(var = "name")

##
print("convert to DESeq dataset")
PURA_ddw <- DESeqDataSetFromSlidingWindows(countData = PURA,
                                           colData = PURA_colData,
                                           annotObj = annotationData,
                                           tidy = FALSE,
                                           design = ~type,
                                           start0based = TRUE)

##
print("filter output based on frequency")
keep <- rowSums(counts(PURA_ddw)) >= 5
PURA_ddw <- PURA_ddw[keep,]
rm(list = "keep")

print("estimate size factors")
PURA_ddw <- estimateSizeFactors(PURA_ddw)

## 
PURA_ddw <- estimateDispersions(PURA_ddw, fitType = "local", quiet = TRUE)
PURA_ddw <- nbinomWaldTest(PURA_ddw)

PURA_resultWindows <- resultsDEWSeq(PURA_ddw,
                                    contrast = c("type", "IP", "SMI"),
                                    parallel = FALSE,
                                    BPPARAM = 12,
                                    tidy = TRUE) %>% as_tibble

PURA_resultWindows[,"p_adj_IHW"] <- adj_pvalues(ihw(pSlidingWindows ~ baseMean, 
                                                    data = PURA_resultWindows,
                                                    alpha = 0.05,
                                                    nfolds = 10))

saveRDS(PURA_ddw, "./DEW-seq/Matrix_w50s20/PURA_ddw_ct5.RDS")
saveRDS(PURA_resultWindows, "./DEW-seq/Matrix_w50s20/PURA_resultWindows.RDS")


pdf("./Plots/DEW-seq_Dispersion_PURA_ct5.pdf", width = 7, height = 7)
plotDispEsts(PURA_ddw)
dev.off()

## extract regions
extract_regions(resultWindows = PURA_resultWindows,
                padjThresh = 0.01,
                log2FoldChangeThresh = 2,
                trackName = "Binding site PURA 0.01",
                fileName = "./DEW-seq/Regions_w50s20/PURA_enrichedregions.log2FC_2.0_p0.01.bed") %>%
  saveRDS("./DEW-seq/Regions_w50s20/PURA_enrichedregions.log2FC_2.0_p0.01.RDS")


extract_regions(resultWindows = PURA_resultWindows,
                padjThresh = 0.05,
                log2FoldChangeThresh = 2,
                trackName = "Binding site PURA 0.05",
                fileName = "./DEW-seq/Regions_w50s20/PURA_enrichedregions.log2FC_2.0_p0.05.bed") %>%
  saveRDS("./DEW-seq/Regions_w50s20/PURA_enrichedregions.log2FC_2.0_p0.05.RDS")



# * * * * * PURB
print("PURB")
##
print("prep data")
PURB <- countData %>%
  dplyr::select(unique_id, contains("PURB"), -contains("GFP"), -contains("PURA")) %>%
  #dplyr::select(unique_id, contains("IP"), contains("SMI")) %>%
  column_to_rownames(var = "unique_id")

PURB_colData <- tibble(name = colnames(PURB)) %>%
  mutate(type = case_when(
    #str_detect(name, "IP") ~ "IP",
    str_detect(name, "SMI") ~ "SMI",
    TRUE ~ "IP"
  )) %>%
  mutate(type = fct_relevel(type, c("IP", "SMI"))) %>%
  column_to_rownames(var = "name")

##
print("convert to DESeq dataset")
PURB_ddw <- DESeqDataSetFromSlidingWindows(countData = PURB,
                                           colData = PURB_colData,
                                           annotObj = annotationData,
                                           tidy = FALSE,
                                           design = ~type,
                                           start0based = TRUE)

##
print("filter output based on frequency")
keep <- rowSums(counts(PURB_ddw)) >= 5
PURB_ddw <- PURB_ddw[keep,]
rm(list = "keep")

print("estimate size factors")
PURB_ddw <- estimateSizeFactors(PURB_ddw)

## 
PURB_ddw <- estimateDispersions(PURB_ddw, fitType = "local", quiet = TRUE)
PURB_ddw <- nbinomWaldTest(PURB_ddw)

PURB_resultWindows <- resultsDEWSeq(PURB_ddw,
                                    contrast = c("type", "IP", "SMI"),
                                    parallel = FALSE,
                                    BPPARAM = 12,
                                    tidy = TRUE) %>% as_tibble

PURB_resultWindows[,"p_adj_IHW"] <- adj_pvalues(ihw(pSlidingWindows ~ baseMean, 
                                                    data = PURB_resultWindows,
                                                    alpha = 0.05,
                                                    nfolds = 10))

saveRDS(PURB_ddw, "./DEW-seq/Matrix_w50s20/PURB_ddw_ct5.RDS")
saveRDS(PURB_resultWindows, "./DEW-seq/Matrix_w50s20/PURB_resultWindows.RDS")


pdf("./Plots/DEW-seq_Dispersion_PURB_ct5.pdf", width = 7, height = 7)
plotDispEsts(PURB_ddw)
dev.off()

## extract regions
extract_regions(resultWindows = PURB_resultWindows,
                padjThresh = 0.01,
                log2FoldChangeThresh = 2,
                trackName = "Binding site PURB 0.01",
                fileName = "./DEW-seq/Regions_w50s20/PURB_enrichedregions.log2FC_2.0_p0.01.bed") %>%
  saveRDS("./DEW-seq/Regions_w50s20/PURB_enrichedregions.log2FC_2.0_p0.01.RDS")


extract_regions(resultWindows = PURB_resultWindows,
                padjThresh = 0.05,
                log2FoldChangeThresh = 2,
                trackName = "Binding site PURB 0.05",
                fileName = "./DEW-seq/Regions_w50s20/PURB_enrichedregions.log2FC_2.0_p0.05.bed") %>%
  saveRDS("./DEW-seq/Regions_w50s20/PURB_enrichedregions.log2FC_2.0_p0.05.RDS")






