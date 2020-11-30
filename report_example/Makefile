DOCUMENT	= report_example
FIGSDIR		= figures
THTMPL          = $(wildcard conf/*.tex)

all: $(DOCUMENT).pdf

#Necessary to force immediate expansion of variable in sub Makefiles
FIGURES		:= $(wildcard figures/*.eps)
INTERMEDIATES	:= .ltarget $(foreach ext,aux bbl blg pdf lof log lot out toc,$(DOCUMENT).$(ext))

MAKEFILES	:= $(wildcard figures/*/Makefile.in)

include $(MAKEFILES)

$(DOCUMENT).pdf: %.pdf: %.bbl
	pdflatex -shell-escape "\def\nogenfigs{}\input{$(<:.bbl=.tex)}"
	pdflatex -shell-escape $(<:.bbl=.tex)

$(DOCUMENT).bbl: .ltarget
	bibtex $(@:.bbl=) || echo "Skipping BibTeX"

.ltarget: $(DOCUMENT).tex $(wildcard *.tex) $(wildcard *.bib) drdc.bst $(FIGURES)
	touch $@ && pdflatex -shell-escape "\def\nogenfigs{}\input{$<}"

.SECONDARY:	$(INTERMEDIATES)

clean:
	rm -f .ltarget *.aux *.bbl *.blg *.log *.dvi *.600pk *.mp *.tfm *.t1 *.out *.toc *.lof *.brf *.lot $(DOCUMENT).ps
	rm -f $(foreach dir,$(dir  $(wildcard figures/*/Makefile.in)),$(dir)*.aux $(dir)*.log $(dir)*.dvi $(dir)*.out $(dir)*.exe $(dir)*.pls)

clear: clean
	rm -r $(DOCUMENT).pdf
	rm -f $(patsubst %/,%.pdf,$(foreach dir,$(dir  $(wildcard figures/*/Makefile.in)),$(dir)))
