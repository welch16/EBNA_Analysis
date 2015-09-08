
clean:
	rm -fr *~
	rm -fr */*~
	rm -f .RDataTmp
	rm -f .Rhistory
	rm -f */.Rhistory
	rm -f .RData
	rm -f .*~
	rm -f *Rout

# knit the vignettes
notes/%.md:rmarkdown/%.Rmd
	cd rmarkdown;R -e 'library(knitr);knit("$(<F)")';mv $(<F:.Rmd=.md) ../notes;cd ..

data/ranges/all_EBV_GenomicRanges.RData:
	R CMD BATCH --no-save scripts/create_GRanges_EBV.R

data/RData/pam_clusters.RData:
	R CMD BATCH --no-save scripts/calculate_pam_cluster.R

data/RData/annotation.RData:
	R CMD BATCH --no-save scripts/annotate_EBV.R

data/RData/annotation_dnase.RData:
	R CMD BATCH --no-save scripts/annotate_EBV_with_Dnase.R

data/RData/association_tests.RData:
	R CMD BATCH --no-save scripts/association_on_annot_tests.R

data/RData/chromHMM_segway_proportions.RData:
	R CMD BATCH --no-save scripts/chrom_hmm_segway_analysis.R

data/RData/chromHMM_segway_proportions_dnase.RData:
	R CMD BATCH --no-save scripts/chrom_hmm_segway_analysis_dnase.R

data/RData/chromHMM_proportions.RData:
	R CMD BATCH --no-save scripts/chrom_hmm_analysis.R

data/RData/chromHMM_proportions_dnase.RData:
	R CMD BATCH --no-save scripts/chrom_hmm_analysis_dnase.R

inst/generated/Big_overlaps_matrix.csv:
	R CMD BATCH --no-save scripts/unify_overlap_matrix.R

data/RData/unified_lists_wProbs.RData:
	/unsup/R-3.2.1/bin/R CMD BATCH --no-save scripts/build_base_lists.R

data/RData/summits.RData:data/RData/unified_lists_wProbs.RData
	/unsup/R-3.2.1/bin/R CMD BATCH --no-save scripts/find_summits.R

data/RData/factor_overlaps.RData:data/RData/unified_lists_wProbs.RData
	/unsup/R-3.2.1/bin/R CMD BATCH --no-save scripts/factor_overlaps.R

data/RData/sequences_around_summit.RData:data/RData/summits.RData data/RData/unified_lists_wProbs.RData
	/unsup/R-3.2.1/bin/R CMD BATCH --no-save scripts/extract_sequences.R

figures/for_paper/fig1.pdf:data/RData/unified_lists_wProbs.RData
	/unsup/R-3.2.1/bin/R CMD BATCH --no-save scripts/fig1.R

figures/for_paper/fig3C_EBNA3B.pdf:data/RData/unified_lists_wProbs.RData data/RData/factor_overlaps.RData
	/unsup/R-3.2.1/bin/R CMD BATCH --no-save scripts/fig3C.R

data/RData/TF_overlap_proportion.RData:data/RData/factor_overlaps.RData data/RData/unified_lists_wProbs.RData
	/unsup/R-3.2.1/bin/R CMD BATCH --no-save scripts/proportion_plots.R

profiles:
	/unsup/R-3.2.1/bin/R CMD BATCH --no-save scripts/histone_profiles.R

data/RData/pam_analysis_K10.RData:
	/unsup/R-3.2.1/bin/R CMD BATCH --no-save scripts/histone_signal_analysis.R

sequences:data/RData/sequences_around_summit.RData data/RData/factor_overlaps.RData
	/unsup/R-3.2.1/bin/R CMD BATCH --no-save scripts/print_sequences.R

figS1:data/RData/factor_overlaps.RData
	/unsup/R-3.2.1/bin/R CMD BATCH --no-save scripts/figS1.R

