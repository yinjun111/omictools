#! /bin/sh

if [ "$USER" !=  "ubuntu" ]; then 
	echo "Script must be run by ubuntu account" 
	exit 1
fi

cd ~/Programs
wget -O - http://cpanmin.us | sudo perl - --self-upgrade
sudo mkdir /apps/perl5lib

# install modules
echo "export PERL5LIB=/apps/perl5lib/lib/perl5" >> /home/ubuntu/.bashrc
export PERL5LIB=/apps/perl5lib/lib/perl5
echo $PERL5LIB

sudo /usr/local/bin/cpanm --local-lib=/apps/perl5lib/ List::Utils
sudo /usr/local/bin/cpanm --local-lib=/apps/perl5lib/ File::Which
sudo /usr/local/bin/cpanm --local-lib=/apps/perl5lib/ Digest::MD5
sudo /usr/local/bin/cpanm --local-lib=/apps/perl5lib/ Excel::Writer::XLSX
sudo /usr/local/bin/cpanm --local-lib=/apps/perl5lib/ Text::CSV

sudo chmod -R 755 /apps/perl5lib/
