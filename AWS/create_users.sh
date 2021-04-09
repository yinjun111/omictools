#!/bin/bash

. "/etc/parallelcluster/cfnconfig"

IFS=","

if [ "${cfn_node_type}" = "ComputeFleet" ]; then
    while read USERNAME USERID
    do
        # -M do not create home since Master node is exporting /homes via NFS
        # -u to set UID to match what is set on the Master node
        if ! [ $(id -u $USERNAME 2>/dev/null || echo -1) -ge 0 ]; then
			#add sudo here ... 
            sudo useradd -M -u $USERID $USERNAME
        fi
    done < "/data/users/userlistfile"
fi
