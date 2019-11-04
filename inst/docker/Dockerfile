FROM pshannon/trenashinybase:latest
USER trena
WORKDIR /home/trena
COPY packages/*gz ./
COPY demoScripts/* ./
RUN R CMD INSTALL -l ~/library TrenaProject_0.99.11.tar.gz
RUN R CMD INSTALL -l ~/library TrenaProjectPlacenta_0.99.1.tar.gz
USER root


