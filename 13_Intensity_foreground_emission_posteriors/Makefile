#DOCNAME=13_Intensity_foreground_emission_posteriors
DOCNAME=BP_XIV_Intensity_foreground_sampling_degeneracies_and_priors_for_AandA

all: paper

.PHONY: clean

paper:
	pdflatex $(DOCNAME).tex
	bibtex $(DOCNAME).aux
	pdflatex $(DOCNAME).tex
	pdflatex $(DOCNAME).tex

view: report
	evince $(DOCNAME).pdf

clean:
	rm *.blg *.bbl *.aux *.log
