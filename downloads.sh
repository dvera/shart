wget -O- https://github.com/lh3/bwa/archive/v0.7.15.tar.gz | tar zx && \
 cd /root/bwa-0.7.15/ && \
 make && \
 cp bwa qualfa2fq.pl xa2multi.pl /usr/local/bin/ && \
 rm -fr /root/bwa-0.7.15

wget -O- https://github.com/arq5x/bedtools2/releases/download/v2.26.0/bedtools-2.26.0.tar.gz | tar zx && \
 cd /root/bedtools2/ && \
 make && \
 cp bin/* /usr/local/bin/ && \
 rm -fr /root/bedtools2/

wget https://github.com/samtools/samtools/releases/download/1.4/samtools-1.4.tar.bz2 && \
 bunzip2 samtools-1.4.tar.bz2 && \
 tar -xvf samtools-1.4.tar && \
 cd /root/samtools-1.4 && \
 ./configure --without-curses && \
 make && \
 make install && \
 rm -f /root/samtools-1.4.tar && \
 rm -fr /root/samtools-1.4

cd /opt && \
 wget http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.5.zip && \
 unzip fastqc_v0.11.5.zip && \
 chmod +x FastQC/fastqc && \
 rm -f fastqc_v0.11.5.zip

pip install cutadapt

R -e 'install.packages("devtools", repos="http://cran.us.r-project.org")' && \
 R -e 'devtools::install_github("dvera/conifur")' && \
 R -e 'devtools::install_github("dvera/converge")' && \
 R -e 'devtools::install_github("dvera/gyro")' && \
 R -e 'devtools::install_github("dvera/travis")'
