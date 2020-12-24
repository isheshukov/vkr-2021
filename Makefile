now      := $(shell date +'%Y-%m-%dT%H-%M')
name-lin := vkr-$(now)
name-win := vkr-$(now)-w64
name-doc := report-$(now)

all: lin win doc

lin:
	git archive --format=zip -o dist/$(name-lin).zip HEAD

win:
	mkdir -p dist/$(name-win)/vkr
	cp -r windata/* dist/$(name-win)/
	git archive HEAD | tar -x -C dist/$(name-win)/vkr
	cd dist/$(name-win) && zip -r -9 ../$(name-win).zip . -x "*share/vim*"
	rm -r dist/$(name-win)

doc:
	cd doc && latexmk -C
	cd doc && latexmk -pdf -shell-escape thesis.tex
	mkdir -p dist/$(name-doc)
	find doc -type f -regex  '.*\(tex\|pdf\|bib\|rtx\)$$' -exec cp '{}' dist/$(name-doc) \;
	cd dist/$(name-doc) && zip -r -9 ../$(name-doc).zip .
	rm -r dist/$(name-doc)

.PHONY: all doc
