all:	pdf clean

pdf:	tr_2005_02.tex
	epstopdf splitSet.eps
	pdflatex tr_2005_02
	pdflatex tr_2005_02
	pdflatex tr_2005_02

clean:
	rm -f *.aux	*.eps *.log *.out *.tex *.toc
	rm -f Rplots.ps splitSet.pdf tr_2005_02-*
