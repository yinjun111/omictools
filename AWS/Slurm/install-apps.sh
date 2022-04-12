#! /bin/sh

if [ "$USER" !=  "ubuntu" ]; then 
	echo "Script must be run by ubuntu account" 
	exit 1
fi

echo "Installing applications"
sudo chmod 755 /apps
sudo chown ubuntu.ubuntu /apps
mkdir /home/ubuntu/Programs
mkdir /home/ubuntu/Pipeline


exit
# Install Zip30
cd ~/Programs
rm -rf ~/Programs/zip30 ~/Programs/zip30.tar.gz
wget https://downloads.sourceforge.net/infozip/zip30.tar.gz
tar xzvf zip30.tar.gz
cd ~/Programs/zip30
make -f unix/Makefile generic
sudo cp -R ~/Programs/zip30 /apps/
