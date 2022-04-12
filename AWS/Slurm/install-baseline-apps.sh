#! /bin/sh

if [ "$USER" !=  "ubuntu" ]; then 
	echo "Script must be run by ubuntu account" 
	exit 1
fi

echo "Updating O/S and installing system software"
sudo apt -y update
sudo apt -y upgrade
sudo apt -y install emacs git wget zip screen
sudo apt -y install gcc gfortran g++ 
sudo apt -y zlib1g-dev zlib1g bzip2 libbz2-dev xz-utils pcre2-utils libpcre2-dev libcurl4-gnutls-dev
sudo apt -y install libcairo2-dev proj-bin libproj-dev
sudo apt -y install libharfbuzz-dev libfribidi-dev
sudo apt -y install libfreetype-dev libpng-dev libtiff-dev libjpeg-turbo8-dev
sudo apt -y install default-jdk

#default parallel in ubuntu doesn't work
sudo apt-get install parallel
