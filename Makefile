all: hw1.pdf

%.pdf: %.tex clean
	pdflatex $<
	pdflatex $<

clean:
	rm *.pdf || true
	rm *.aux || true
	rm *.out || true

.PHONY: clean
