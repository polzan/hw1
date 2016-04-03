all: hw1.pdf hw1_code.pdf

hw1.pdf: hw1_code.pdf

%.pdf: %.tex clean
	pdflatex $<
	pdflatex $<

clean:
	rm *.pdf || true
	rm *.aux || true
	rm *.out || true

.PHONY: clean all
