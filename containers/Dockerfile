FROM quay.io/biocontainers/bamutil:1.0.14--hba45d7a_4

COPY seqkit_linux_amd64.tar.gz /tmp

RUN cd /tmp \
	&& gunzip seqkit_linux_amd64.tar.gz \
	&& tar xf seqkit_linux_amd64.tar \
	&& chmod +x seqkit && mv seqkit /usr/bin/

