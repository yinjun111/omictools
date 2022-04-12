#! /bin/bash
#! /bin/sh

if [ "$USER" !=  "ubuntu" ]; then 
	echo "Script must be run by ubuntu account" 
	exit 1
fi

cd ~/Pipeline/
/usr/local/bin/aws s3 sync s3://ferring-omictools/omictools omictools/
chmod 755 -R ~/Pipeline/omictools/
cd omictools
sh cp_scripts.sh
sudo ln -s /apps/omictools/omictools_caller.pl /usr/bin/omictools

cd ~/Pipeline/
/usr/local/bin/aws s3 sync s3://ferring-omictools/gstools gstools/
chmod 755 -R ~/Pipeline/gstools/
cd gstools/
sh cp_scripts_ubuntu.sh
sudo ln -s /apps/gstools/gstools_caller.pl /usr/bin/gstools

# copy test data
mkdir -p /data/jyin/Pipeline_Test
sudo chown -R ubuntu.ubuntu /data/jyin/Pipeline_Test
/usr/local/bin/aws s3 sync s3://ferring-omictools/Testdata/ /data/jyin/Pipeline_Test/Testdata/

# copy GATK3 data
aws s3 cp s3://ferring-omictools/Others/GATK3/ /apps/GATK3/ --recursive

#copy database
aws s3 sync s3://ferring-omictools/Databases/ /data/jyin/Databases/ --delete

