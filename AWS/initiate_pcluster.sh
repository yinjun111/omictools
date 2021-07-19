#!/bin/sh

#may not work

name=$1
config=$2
key=$3

#create pcluster
#pcluster create $name -c $config

echo "pcluster $name created!"
date

#scp the initiate script
pclusterpublicip=$(pcluster status $name | grep MasterPublicIP | awk '{print $2}')

echo "pcluster $name public IP: $pclusterpublicip"

#add key
#ssh-keyscan $pclusterpublicip >> /home/centos/.ssh/known_hosts


#doesn't work here, because iam doesn't have s3 permission

#run initiate_master node
ssh -t centos@$pclusterpublicip -i $key  << EOF
	hostname
	#aws s3 cp s3://ferring-omictools/omictools/AWS/initiate_masternode.sh /home/centos/
	sudo sh /home/centos/initiate_masternode.sh
	sudo cat /etc/passwd
	
EOF

echo "pcluster $name initiated for omictools!"
date
