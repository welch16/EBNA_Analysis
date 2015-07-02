
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

data/RData/association_tests.RData:
	R CMD BATCH scripts/association_on_annot_tests.R
