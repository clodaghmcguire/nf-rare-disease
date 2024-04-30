# Base Image
FROM rstudio/r-base:4.2.1-focal

# Metadata
LABEL base.image="r-base:4.2.1"
LABEL software="ExomeDepth"
LABEL description="CNV calling for DNA capture assays"
LABEL website="https://CRAN.R-project.org/package=ExomeDepth"
LABEL documentation="https://CRAN.R-project.org/package=ExomeDepth"
LABEL license="https://CRAN.R-project.org/package=ExomeDepth"

# Maintainer
LABEL maintainer = KCH Bioinformatics <kch-tr.KCHBioinformatics@nhs.net>

# set working directory
WORKDIR /root
ENV HOME /root

RUN apt-get -y update
RUN apt-get --yes install \
	libcurl4-openssl-dev \
	time \
	libpoppler-cpp-dev \
	texlive \
	texlive-latex-extra \
	texlive-fonts-extra \
	libssl-dev \
	libfontconfig1-dev \
	libxml2-dev

# LateX packages (Locked on 2020 release)
RUN tlmgr init-usertree
RUN tlmgr install framed && \
	tlmgr install mdframed && \
	tlmgr install zref && \
	tlmgr install needspace && \
	tlmgr install dingbat && \
	tlmgr install booktabs && \
	tlmgr install makecell && \
	tlmgr install colortbl && \
	tlmgr install parskip && \
	tlmgr install fancyhdr && \
	tlmgr install lastpage && \
	tlmgr install xcolor

# Set CRAN mirror
RUN echo 'options(repos=structure(c(CRAN="http://cran.ma.imperial.ac.uk/")))' > /root/.Rprofile

# R packages (BioC)
RUN Rscript -e "install.packages('BiocManager');BiocManager::install(version = '3.15', ask = FALSE)" && \
	Rscript -e "BiocManager::install(c('Biostrings','IRanges','Rsamtools','GenomicRanges','GenomicAlignments'), ask = FALSE, update = FALSE)"

# R packages (CRAN)
RUN Rscript -e "install.packages('optparse')" && \
    Rscript -e "install.packages('randomForest')" && \
    Rscript -e "install.packages('aod')" && \
    Rscript -e "install.packages('ExomeDepth')" && \
    Rscript -e "install.packages('plyr')" && \
    Rscript -e "install.packages('knitr')" && \
    Rscript -e "install.packages('kableExtra')" && \
    Rscript -e "install.packages('tinytex')" && \
    Rscript -e "install.packages('xtable')"

# Add exomedepth scripts and report template
ADD ed2vcf.R /root
ADD readCount.R /root
ADD exomeDepth.R /root
ADD exomeDepth.Rnw /root

# Add QC configs
COPY qc_config/* /root/

CMD ["/bin/bash"]