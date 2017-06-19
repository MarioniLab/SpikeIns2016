all: Supplemental_Materials.pdf prickle.pdf

prickle.pdf: prickle.tex refnorm.bib \
	pics/plate_setup.pdf ../real/pics/variance_exp.pdf ../simulations/pics/setplot.pdf ../simulations/clustering/clust_effect.pdf
	pdflatex prickle.tex
	biber prickle
	pdflatex prickle.tex
	pdflatex prickle.tex

Supplemental_Materials.pdf: supplement.tex refnorm.bib \
	../real/pics/qq_separate.pdf ../real/pics/qq_premixed.pdf ../real/pics/total_ercc.pdf ../real/pics/total_sirv.pdf \
	../real/pics/variance_order.pdf ../sequence_check/biophysical/comparison.pdf ../simulations/depth/depth_effect.pdf 
	pdflatex supplement.tex
	biber supplement
	pdflatex supplement.tex
	pdflatex supplement.tex
	mv supplement.pdf $@

