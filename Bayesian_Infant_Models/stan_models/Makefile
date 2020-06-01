all:
	R -e "rmarkdown::render('stan_plot_models.Rmd')"

continuous:
	while :; do inotifywait -e modify *.Rmd; make; done

.PHONY: continuous
