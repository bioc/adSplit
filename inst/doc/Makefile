all:	pdf clean

pdf:	tr_2005_02.tex
	epstopdf splitSet.eps
	pdflatex tr_2005_02
	pdflatex tr_2005_02
	pdflatex tr_2005_02

clean:
	/bin/rm -f *.out *.bbl *.log *.aux *.blg *.toc  *~ err.rda splitSet*
	/bin/rm -f *.dvi *.ps tr_2005_02-* *.html Rplots* tr_2005_02.tex

