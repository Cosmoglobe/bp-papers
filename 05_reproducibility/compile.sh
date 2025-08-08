#!/bin/bash 

pdflatex main.tex 
bibtex main.aux 
pdflatex main.tex 
pdflatex main.tex 

rm BP05_authors_openjournal.aux main.aux main.bbl main.blg main.log main.out mainNotes.bib

open main.pdf
