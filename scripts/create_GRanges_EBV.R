
##################################################################################
#
#  Script_GRanges_EBNA.R	
# 
##################################################################################

rm(list = ls())
library(GenomicRanges)

# Convert our final lists into GRanges objects with summits 

# Final lists on /scratch/ReneWelch/EBNA/Generated/FinalLists
#	 EBNA2_peaks_wfilter_woverlap_wannot_wencode_wEncodeOverlap_June2012_Grouped.txt
#	 EBNA3A_peaks_wfilter_woverlap_wannot_wencode_wEncodeOverlap_June2012_Grouped.txt
#	 EBNA3B_peaks_wfilter_woverlap_wannot_wencode_wEncodeOverlap_June2012_Grouped.txt
#	 EBNA3C_peaks_wfilter_woverlap_wannot_wencode_wEncodeOverlap_June2012_Grouped.txt
#	 JK234_peaks_wfilter_woverlap_wannot_wencode_wEncodeOverlap_June2012_Grouped.txt
#	 JK92_peaks_wfilter_woverlap_wannot_wencode_wEncodeOverlap_June2012_Grouped.txt

# Peaks on /scratch/ReneWelch/EBNA/Generated/peakRanges
#	EBNA2_peaks.RData
#	EBNA3A_peaks.RData
#	EBNA3B_peaks.RData
#	EBNA3C_peaks.RData
#	JK234_peaks.RData
#	JK92_peaks.RData

sets = c("EBNA2","EBNA3A","EBNA3B","EBNA3C","JK92","JK234")
histones = c("H3k27me3","H3k36me3","H3k4me3","H3k9ac","H3k9me3")

# We are on R directory of EBNA folder
listDir <- "data/lists"
rangeDir <- "data/RData"

directories <- c(l = listDir, r = rangeDir)

#-------------------------------------------------------------------------------
build.GRanges <- function(peak, directories)
{
    # Load the infor in the directories and build a GRanges object
    peakData = read.table(file = file.path(directories[1],paste0(peak,"_peaks_wfilter_woverlap_wannot_wencode_wEncodeOverlap_June2012_Grouped.txt")),
      header= TRUE)
    load(file = file.path(directories[2],paste0(peak,"_peaks.RData")))
    peakRange = peaklist
    drivers = names(peakData[29:dim(peakData)[2]])
    cols = c("region","GeneID","distance",names(peakData)[12:17],drivers)
    aux = lapply( cols,function(x,peakData) peakData[ , x],peakData )
    cols[5:6] = c("JK92","JK234")  # Correction to make consistent set names
    names(aux) = cols
    for(k in cols){
      peakRange@elementMetadata@listData[[k]] <- as.matrix(aux[[k]])
    }
    # Add two additional
    peakRange$RBPJ = ifelse(peakRange$JK92 == 1 | peakRange$JK234 == 1,1,0)
    allTF = drivers[ !grepl("H2Z",drivers)  & !grepl("H3K",drivers) & !grepl("H4K",drivers) & !grepl("Dnase",drivers) ]
    peakRange$allTF = ifelse(rowSums(as.matrix(peakData[,allTF])) > 0 ,1,0)
    return(peakRange)
}
#--------------------------------------------------------
filter.outlier <- function(histone,peak,idx,all.profiles)
{
    profiles = all.profiles[[peak]][[histone]]
    profiles = profiles[,as.logical(idx)]
    return(profiles)
}

# Build the GRanges sets
ranges <- lapply(sets,FUN = build.GRanges, directories)
names(ranges) = sets

# We are using the rule:
# If a peak is common to JK92 and JK234 use the JK234 version

jk92 = ranges[["JK92"]]
jk234 = ranges[["JK234"]]
jk92_sub = subset(jk92, JK234 == 0)
ranges[["RBPJ"]] = c (jk234,jk92_sub)

## all.profiles = load.Profiles(sets,histones,profileDir)

# Remove the outlier peak for H3k27me3 and EBNA3C
# The info of the peak is: chr2 33141400 - 33141799
idx = !(seqnames(ranges[["EBNA3C"]]) == "chr2" & start(ranges[["EBNA3C"]])==33141400)
ranges[["EBNA3C"]] = subset(ranges[["EBNA3C"]],idx)

save(file = file.path("data/ranges","all_EBV_GenomicRanges.RData"), list = "ranges")

