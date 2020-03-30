all:
	R -e "rmarkdown::render('plot.Rmd')"

continuous:
	while :; do inotifywait -e modify *.Rmd; make; done

.PHONY: continuous
