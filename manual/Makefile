# Build of the manual.  make pdf     : the comparison table is regenerated from ex/r/*.csv, the PDF is built (3 passes)
#                      make examples: reruns the do-files of ex/ in Stata 17 (needs NWCORE=/path/to/nwcommands; about 5 minutes; the data come from the web)
#                      make compare : the comparison with goldfish of chapter 11 (needs Stata, and R with goldfish 1.6.12)
#                      make check   : undefined references, overfull boxes, glued words
# The PDF is committed with git add -f (the LaTeX build files are ignored).
all: pdf

tables:
	python3 build_compare.py

pdf: tables
	pdflatex -interaction=nonstopmode nwdynam_manual.tex
	pdflatex -interaction=nonstopmode nwdynam_manual.tex
	pdflatex -interaction=nonstopmode nwdynam_manual.tex

# the numerical comparison with goldfish of chapter 11: Stata writes the two event logs, R fits the models, Stata fits them, python builds the table
compare:
	NWCORE=$(NWCORE) ./run_examples.sh ex90_export_for_r
	cd ex/r && Rscript compare.R
	NWCORE=$(NWCORE) ./run_examples.sh ex91_compare_stata
	python3 build_compare.py

examples:
	./run_examples.sh ex01_first ex02_effects ex03_special ex04_post ex05_submodels ex06_opp ex07_limits

check:
	python3 check_layout.py

clean:
	rm -f *.aux *.log *.toc *.out *.bbl *.blg *.lof *.lot
