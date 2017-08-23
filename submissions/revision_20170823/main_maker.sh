cat ../../manuscript/prickle.tex | sed "s/pics\/plate_setup.pdf/Figure1.pdf/" |
    sed "s/..\/real\/pics\/variance_exp.pdf/Figure2.pdf/" |
    sed "s/..\/simulations\/pics\/setplot.pdf/Figure3.pdf/" |
    sed "s/..\/simulations\/clustering\/clust_effect.pdf/Figure4.pdf/" > LunEtAl_manuscript.tex

cp ../../manuscript/pics/plate_setup.pdf Figure1.pdf
cp ../../manuscript/../real/pics/variance_exp.pdf Figure2.pdf
cp ../../manuscript/../simulations/pics/setplot.pdf Figure3.pdf
cp ../../manuscript/../simulations/clustering/clust_effect.pdf Figure4.pdf

cp ../../manuscript/refnorm.bib .

pdflatex LunEtAl_manuscript.tex 
biber LunEtAl_manuscript
pdflatex LunEtAl_manuscript.tex 
pdflatex LunEtAl_manuscript.tex 

