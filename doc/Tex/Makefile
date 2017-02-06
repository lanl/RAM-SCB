# This make file generates the documentation for RAM_SCB
include ../../Makefile.def
include ${DIR}/share/build/Makefile.doc

help:
	@echo "       Make manuals for RAM_SCB"
	@echo ""
	@echo "Targets:"
	@echo ""
	@echo "PDF          - produces PDF manuals in the ../ directory"
	@echo ""
	@echo "clean        - remove temporary files"
	@echo ""
	@echo "cleanpdf     - clean + remove PDF manuals"


RAM_SCBXML:
	@rm -f RAM_SCBXML.tex
	@touch RAM_SCBXML.tex
	@perl -e 'for $$i (glob("../../PARAM.XML ../../[A-Z][A-Z]/*/PARAM.XML")){`${SCRIPTDIR}/XmlToTex.pl $$i >> RAM_SCBXML.tex`}'
	@echo "${SCRIPTDIR}/XmlToTex.pl */PARAM.XML > RAM_SCBXML.tex"

RAM_SCB.dvi:        \
	RAM_SCBXML  \
	RAM_SCB.tex \
	HEADER.tex  \
	RAM_SCB_intro.tex      \
	RAM_SCB_quickstart.tex \
	RAM_SCB_install.tex    \
	RAM_SCB_swmf.tex       \
	RAM_SCB_run.tex


PDF_FILES = RAM_SCB.pdf
PDF: ${PDF_FILES}

clean:
	rm -f *.out *.aux *.log *.toc *.ind *.idx *.ilg *.dvi *.ps *~ *.protex
	rm -f RAM_SCBXML.tex 
	@echo 'rm -f RAM_SCBREF RAM_SCBXML'
	@rm -f ${RAM_SCBXML}

distclean: clean
	cd ..; rm -f ${PDF_FILES}
