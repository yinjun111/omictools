#! /bin/sh

if [ "$USER" !=  "ubuntu" ]; then 
	echo "Script must be run by ubuntu account" 
	exit 1
fi

cd ~/Programs
wget https://cran.r-project.org/src/base/R-4/R-4.0.2.tar.gz
tar xzvf R-4.0.2.tar.gz
cd R-4.0.2
./configure --prefix=/apps/R-4.0.2 --with-x=no
make
sudo make install

#locfit lib
cd ~/Programs 
wget https://cran.r-project.org/src/contrib/Archive/locfit/locfit_1.5-9.4.tar.gz

sudo /apps/R-4.0.2/bin/Rscript ~/Pipeline/omictools/AWS/initiate_rpkgs.R


wget https://repo.anaconda.com/archive/Anaconda3-2020.11-Linux-x86_64.sh
sudo sh Anaconda3-2020.11-Linux-x86_64.sh -b -p /apps/anaconda3

# install Anaconda packages
sudo /apps/anaconda3/bin/python3 -m pip install cutadapt
sudo /apps/anaconda3/bin/python3 -m pip install multiqc
sudo /apps/anaconda3/bin/python3 -m pip install deeptools
sudo /apps/anaconda3/bin/python3 -m pip install RSeQC
