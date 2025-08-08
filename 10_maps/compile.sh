#!/bin/bash 

pdflatex BP_maps.tex && bibtex BP_maps.aux && pdflatex BP_maps.tex && open BP_maps.pdf
