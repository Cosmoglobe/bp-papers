#!/bin/bash 

pdflatex BP05_reproducibility.tex && bibtex BP05_reproducibility.aux && pdflatex BP05_reproducibility.tex && open BP05_reproducibility.pdf
