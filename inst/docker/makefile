build:
	docker build -t trena_placenta ./

bash:
	docker run -ti --rm trena_placenta bash

host.install:
	R -f installTrenaPrequisites.R

rstudio.bioc:
	docker run -e PASSWORD=password -p 8787:8787 bioconductor/release_core2

rstudio.trena:
	docker run -e PASSWORD=password -p 8787:8787 trena_placenta


