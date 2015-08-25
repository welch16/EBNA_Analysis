
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
	R CMD BATCH scripts/create_GRanges_EBV.R

data/RData/pam_clusters.RData:
	R CMD BATCH scripts/calculate_pam_cluster.R

data/RData/annotation.RData:
	R CMD BATCH scripts/annotate_EBV.R

data/RData/annotation_dnase.RData:
	R CMD BATCH scripts/annotate_EBV_with_Dnase.R

data/RData/association_tests.RData:
	R CMD BATCH scripts/association_on_annot_tests.R

data/RData/chromHMM_segway_proportions.RData:
	R CMD BATCH scripts/chrom_hmm_segway_analysis.R

data/RData/chromHMM_segway_proportions_dnase.RData:
	R CMD BATCH scripts/chrom_hmm_segway_analysis_dnase.R

data/RData/chromHMM_proportions.RData:
	R CMD BATCH scripts/chrom_hmm_analysis.R

data/RData/chromHMM_proportions_dnase.RData:
	R CMD BATCH scripts/chrom_hmm_analysis_dnase.R

inst/generated/Big_overlaps_matrix.csv:
	R CMD BATCH scripts/unify_overlap_matrix.R


data/RData/unified_lists_wProbs.RData:
	R CMD BATCH scripts/build_base_lists.R

data/RData/summits.RData:data/RData/unified_lists_wProbs.RData
	R CMD BATCH scripts/find_summits.R

figures/for_paper/fig1.pdf:data/RData/unified_lists_wProbs.RData
	R CMD BATCH scripts/fig1.R

