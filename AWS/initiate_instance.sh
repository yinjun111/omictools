
#script to creat pcluster

#update yum
sudo yum -y update

sudo yum -y install emacs

#install python
mkdir ~/Programs
cd ~/Programs
wget https://repo.anaconda.com/archive/Anaconda3-2020.11-Linux-x86_64.sh

sudo mkdir /apps/
sudo sh Anaconda3-2020.11-Linux-x86_64.sh -b -p /apps/anaconda3

#install parallelcluster
sudo /apps/anaconda3/bin/pip3 install aws-parallelcluster --upgrade

#configure aws
aws configure

#start cluster #cluster can't have the same name with previous ones
pcluster create torque01
