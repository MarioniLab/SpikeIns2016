cat ../../manuscript/prickle.tex | sed -e '1,/\\end{quote}/d' | cat preamble.tex - | \
     sed "s/section{/section*{/g" | sed "s/bibliographystyle{plainnat}/bibliographystyle{abbrv}/" | 
     sed -E "s/(includegraphics.*\{)/\1..\/..\/manuscript\//" | sed "s/bibliography{refnorm}/bibliography{..\/..\/manuscript\/refnorm}/" |
     sed "s/\\\cite{risso2014normalization}/Risso \\\textit{et al.} \\\cite{risso2014normalization}/" |
     sed "s/\\\cite{brennecke2013accounting}/Brennecke \\\textit{et al.} \\\cite{brennecke2013accounting}/" |
     sed "s/citep{/cite{/g" \
     > final.tex
pdflatex final.tex
bibtex final
pdflatex final.tex
pdflatex final.tex

cp ../../manuscript/supplement.pdf supp.pdf
