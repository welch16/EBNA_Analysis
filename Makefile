
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
