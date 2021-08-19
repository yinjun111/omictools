#update yum

echo "Initiating master node for omictools.\n"
date

sudo yum -y update

sudo yum -y install emacs git wget zip screen

#####
#Directories for shared apps and data
#####

sudo chmod 755 /apps
sudo chown centos /apps

mkdir /home/centos/Programs
sudo chown centos /home/centos/Programs/

mkdir /home/centos/Pipeline/
sudo chown centos /home/centos/Pipeline/

#create users
mkdir /data/users

for user in jyin qmeng yyue mparks yxu cbadger;do
	sudo useradd $user
	echo "$user,`id -u $user`" >> /data/users/userlistfile
done


#####
#System tools
#####

#zip
cd /home/centos/Programs/
wget https://downloads.sourceforge.net/infozip/zip30.tar.gz
tar xzvf zip30.tar.gz

cd /home/centos/Programs/zip30
make -f unix/Makefile generic
sudo cp -R /home/centos/Programs/zip30 /apps/

#parallel
cd /home/centos/Programs/
wget http://ftp.gnu.org/gnu/parallel/parallel-20210322.tar.bz2
tar xjvf parallel-20210322.tar.bz2
cd /home/centos/Programs/parallel-20210322
./configure && make && sudo make install


########
#Omictools
########

#copy omicools
cd /home/centos/Pipeline/
/usr/local/bin/aws s3 sync s3://ferring-omictools/omictools omictools/
chmod 755 -R /home/centos/Pipeline/omictools/

cd omictools/
sh cp_scripts.sh
sudo ln -s /apps/omictools/omictools_caller.pl /usr/bin/omictools

#copy gstool
cd /home/centos/Pipeline/
/usr/local/bin/aws s3 sync s3://ferring-omictools/gstools gstools/
chmod 755 -R /home/centos/Pipeline/gstools/

cd gstools/
sh cp_scripts.sh
sudo ln -s /apps/gstools/gstools_caller.pl /usr/bin/gstools

#copy database
/usr/local/bin/aws s3 sync s3://ferring-omictools/Databases/ /data/jyin/Databases/ --delete

#copy testdata
mkdir -p /data/jyin/Pipeline_Test
sudo chown -R centos /data/jyin/Pipeline_Test
/usr/local/bin/aws s3 sync s3://ferring-omictools/Testdata/ /data/jyin/Pipeline_Test/Testdata/

#copy GATK3
aws s3 cp s3://ferring-omictools/Others/GATK3/ /apps/GATK3/ --recursive

#####
#install R #need to remove anaconda from PATH
#####

#https://cran.r-project.org/src/base/R-4/R-4.0.4.tar.gz
#https://cran.r-project.org/src/base/R-4/R-4.0.2.tar.gz


cd /home/centos/Programs

sudo yum -y install gcc gcc-gfortran gcc-c++ zlib-devel bzip2 bzip2-devel xz-devel pcre2 pcre2-devel curl-devel

wget https://cran.r-project.org/src/base/R-4/R-4.0.2.tar.gz
tar xzvf R-4.0.2.tar.gz
cd /home/centos/Programs/R-4.0.2
./configure --prefix=/apps/R-4.0.2 --with-x=no
make
sudo make install

#for Cairo
sudo yum -y install cairo cairo-devel proj proj-devel
sudo yum -y install harfbuzz-devel fribidi-devel
sudo yum -y install freetype-devel libpng-devel libtiff-devel libjpeg-turbo-devel

#Install R packages
sudo /apps/R-4.0.2/bin/Rscript /home/centos/Pipeline/omictools/AWS/initiate_rpkgs.R

#########
#download anaconda, then choose /apps/anaconda3
##########
cd /home/centos/Programs
wget https://repo.anaconda.com/archive/Anaconda3-2020.11-Linux-x86_64.sh
sudo sh Anaconda3-2020.11-Linux-x86_64.sh -b -p /apps/anaconda3

#install cutadapt
sudo /apps/anaconda3/bin/python3 -m pip install cutadapt

#install multiqc
sudo /apps/anaconda3/bin/python3 -m pip install multiqc

#Bamcoverage
sudo /apps/anaconda3/bin/python3 -m pip install deeptools

#RSeQC
sudo /apps/anaconda3/bin/python3 -m pip install RSeQC

#######
#Install perl modules
#######

sudo yum -y install cpan

#install cpanm
wget -O - http://cpanmin.us | sudo perl - --self-upgrade
 
#install modules
mkdir /apps/perl5lib/

echo "export PERL5LIB=/apps/perl5lib/lib/perl5" >> /home/centos/.bashrc;source /home/centos/.bashrc

/usr/local/bin/cpanm --local-lib=/apps/perl5lib/ List::Utils
/usr/local/bin/cpanm --local-lib=/apps/perl5lib/ File::Which
/usr/local/bin/cpanm --local-lib=/apps/perl5lib/ Digest::MD5
/usr/local/bin/cpanm --local-lib=/apps/perl5lib/ Excel::Writer::XLSX
/usr/local/bin/cpanm --local-lib=/apps/perl5lib/ Text::CSV

sudo chmod -R 755 /apps/perl5lib/

########
#Programs for RNASeq
########

cd /home/centos/Programs
sudo chown centos /home/centos/Pipeline/

#Fastqc
wget https://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.9.zip
unzip fastqc_v0.11.9.zip
sudo cp -R FastQC /apps/
chmod +x /apps/FastQC/fastqc

#RSEM
wget https://github.com/deweylab/RSEM/archive/v1.3.3.tar.gz
tar xzvf v1.3.3.tar.gz
cd RSEM-1.3.3/
make
sudo cp -R /home/centos/Programs/RSEM-1.3.3/ /apps/ 
 
#STAR
cd /home/centos/Programs
wget https://github.com/alexdobin/STAR/archive/2.7.8a.tar.gz
tar xzvf 2.7.8a.tar.gz
sudo cp -R STAR-2.7.8a/ /apps/
#need to recompile?
cd /apps/STAR-2.7.8a/source
sudo make STAR

#install samtools
cd /home/centos/Programs
wget https://github.com/samtools/samtools/releases/download/1.12/samtools-1.12.tar.bz2
tar xjvf samtools-1.12.tar.bz2
cd /home/centos/Programs/samtools-1.12
./configure --prefix=/apps/samtools-1.12
make
sudo make install

#Homer
mkdir /apps/Homer
cd /apps/Homer
wget http://homer.ucsd.edu/homer/configureHomer.pl
perl configureHomer.pl -install
 
#GSEA
cd /home/centos/Programs
wget https://data.broadinstitute.org/gsea-msigdb/gsea/software/desktop/4.0/GSEA_Linux_4.0.3.zip
unzip GSEA_Linux_4.0.3.zip
sudo cp -R GSEA_Linux_4.0.3 /apps/

#sratools
cd /home/centos/Programs
wget  https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/2.11.0/sratoolkit.2.11.0-centos_linux64.tar.gz

tar xzvf sratoolkit.2.11.0-centos_linux64.tar.gz

sudo cp -R sratoolkit.2.11.0-centos_linux64/ /apps/

#htslib for GATK3
cd /home/centos/Programs
wget https://github.com/samtools/htslib/releases/download/1.13/htslib-1.13.tar.bz2

cd /home/centos/Programs/htslib-1.13
./configure --prefix=/apps/htslib-1.13
make
sudo make install

#snpEff
cd /home/centos/Programs
wget https://snpeff.blob.core.windows.net/versions/snpEff_latest_core.zip

unzip snpEff_latest_core.zip
cd /home/centos/Programs/snpEff
sudo cp -R /home/centos/Programs/snpEff /apps/snpEff

sudo chmod 777 /apps/snpEff

java -jar /apps/snpEff/snpEff.jar download hg38
java -jar /apps/snpEff/snpEff.jar download mm10

#######
#chown for some folders
#######

sudo chown -R centos /data/jyin/
sudo chown -R centos /home/centos/Programs
sudo chown -R centos /home/centos/Pipeline
sudo chown -R centos /apps/perl5lib/

echo "Pcluster initiation completed!"
date

