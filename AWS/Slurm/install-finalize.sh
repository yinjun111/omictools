#! /bin/sh

if [ "$USER" !=  "ubuntu" ]; then 
	echo "Script must be run by ubuntu account" 
	exit 1
fi

sudo chown -R ubuntu.ubuntu /apps/* 
sudo chown -R ubuntu.ubuntu /data/jyin

