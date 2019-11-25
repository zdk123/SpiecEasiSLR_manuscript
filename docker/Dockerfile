## Docker container recipe for reproducing results in this repo
FROM continuumio/miniconda3:4.7.12
LABEL maintainer="Zachary D Kurtz <zdkurtz@gmail.com>"

COPY environment.yml ./
RUN conda env update -n base --file environment.yml --prune

RUN ln -s /bin/tar /bin/gtar && \
    Rscript -e "install.packages('https://cran.r-project.org/src/contrib/stringi_1.4.3.tar.gz')" && \
    Rscript -e "devtools::install_github('zdk123/SpiecEasi', ref='d6cf2fb', upgrade='never')"

ENTRYPOINT [ "/opt/conda/bin/tini", "--" ]
CMD [ "/bin/bash" ]
