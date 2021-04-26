now      := $(shell date +'%Y-%m-%dT%H-%M')
name-lin := vkr-$(now)
name-win := vkr-$(now)-w64
name-thesis := thesis-$(now)
name-slides := slides-$(now)

all: lin win doc
doc: thesis slides

lin:
	git archive --format=zip -o dist/$(name-lin).zip HEAD

win:
	mkdir -p dist/$(name-win)/vkr
	cp -r windata/* dist/$(name-win)/
	git archive HEAD | tar -x -C dist/$(name-win)/vkr
	cd dist/$(name-win) && zip -r -9 ../$(name-win).zip . -x "*share/vim*"
	rm -r dist/$(name-win)

thesis:
	cd doc/thesis && latexmk -C
	cd doc/thesis && latexmk -pdf -shell-escape thesis.tex
	mkdir -p dist/$(name-thesis)
	find doc/thesis -type f -regex  '.*\(tex\|pdf\|bib\|rtx\)$$' -exec cp '{}' dist/$(name-thesis) \;
	cd dist/$(name-thesis) && zip -r -9 ../$(name-thesis).zip .
	rm -r dist/$(name-thesis)
slides:
	cd doc/pres && latexmk -C
	cd doc/pres && latexmk -pdf -shell-escape main.tex
	mkdir -p dist/$(name-slides)
	find doc/pres -type f -regex  '.*\(tex\|pdf\|bib\|rtx\)$$' -exec cp '{}' dist/$(name-slides) \;
	cd dist/$(name-slides) && zip -r -9 ../$(name-slides).zip .
	rm -r dist/$(name-slides)

.PHONY: all
