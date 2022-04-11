#synchronize scripts in the dev codes folder with the running folder
rsync -a -d --delete /home/jyin/Pipeline/omictools/ /apps/omictools/

/usr/local/bin/aws s3 sync /home/jyin/Pipeline/omictools/ s3://ferring-omictools/omictools --delete

#make everything execuable +x
chmod 755 -R /apps/omictools/


