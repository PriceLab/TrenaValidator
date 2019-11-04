all:  docs install

docs:
	R -e "devtools::document()"
build:
	(cd ..; R CMD build --no-build-vignettes trenaValidator)

install:
	(cd ..; R CMD INSTALL --no-test-load trenaValidator)

check:
	(cd ..; R CMD check `ls -t trenaValidator) | head -1`)

test:
	for x in inst/unitTests/test_*.R; do echo $$x; R -f $$x; done

