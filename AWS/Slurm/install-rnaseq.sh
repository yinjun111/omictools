#! /bin/sh

if [ "$USER" !=  "ubuntu" ]; then 
	echo "Script must be run by ubuntu account" 
	exit 1
fi

# Fastqc
cd ~/Programs
wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip
unzip fastqc_v0.11.9.zip
sudo cp -R FastQC /apps/
sudo chmod +x /apps/FastQC/fastqc

#RSEM
cd ~/Programs
wget https://github.com/deweylab/RSEM/archive/v1.3.3.tar.gz
tar xzvf v1.3.3.tar.gz
cd RSEM-1.3.3/
make
sudo cp -R ~/Programs/RSEM-1.3.3/ /apps/ 

# STAR
cd ~/Programs
wget https://github.com/alexdobin/STAR/archive/2.7.8a.tar.gz
tar xzvf 2.7.8a.tar.gz
sudo cp -R STAR-2.7.8a/ /apps/
#need to recompile?
cd /apps/STAR-2.7.8a/source
sudo make STAR
sudo cp STAR ../bin/Linux_x86_64

# samtools
cd ~/Programs
wget https://github.com/samtools/samtools/releases/download/1.12/samtools-1.12.tar.bz2
tar xjvf samtools-1.12.tar.bz2
cd samtools-1.12
./configure --prefix=/apps/samtools-1.12
make
sudo make install

#Homer
sudo mkdir /apps/Homer
cd /apps/Homer
sudo wget http://homer.ucsd.edu/homer/configureHomer.pl
sudo perl configureHomer.pl -install

#GSEA
cd ~/Programs
# note: need to download the software from http://www.gsea-msigdb.org after creating an account
unzip GSEA_Linux_4.2.3.zip
sudo cp -R GSEA_Linux_4.2.3 /apps/

# sratools
cd ~/Programs
wget  https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.11.0/sratoolkit.2.11.0-centos_linux64.tar.gz
tar xzvf sratoolkit.2.11.0-centos_linux64.tar.gz
sudo cp -R sratoolkit.2.11.0-centos_linux64/ /apps/

#htslib for GATK3
cd ~/Programs
wget https://github.com/samtools/htslib/releases/download/1.13/htslib-1.13.tar.bz2
tar xf htslib-1.13.tar.bz2
cd htslib-1.13
./configure --prefix=/apps/htslib-1.13
make
sudo make install

# snpEff
cd ~/Programs
cd ~/Programs
wget https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip

unzip snpEff_latest_core.zip
sudo cp -R snpEff /apps/snpEff

sudo chmod 777 /apps/snpEff

java -jar /apps/snpEff/snpEff.jar download hg38
java -jar /apps/snpEff/snpEff.jar download mm10

#featureCounts
cd ~/Programs
wget https://downloads.sourceforge.net/project/subread/subread-2.0.3/subread-2.0.3-Linux-x86_64.tar.gz --no-check-certificate
tar xzvf subread-2.0.3-Linux-x86_64.tar.gz
sudo cp -R subread-2.0.3-Linux-x86_64 /apps/

