#! /bin/bash

for user in jyin qmeng yxu cbadger jpatterson; 
do 
   # sudo adduser $user
   # sudo passwd -e $user
   #sudo mkdir /home/$user/.ssh
   echo $user
done

sudo echo "ssh-rsa AAAAB3NzaC1yc2EAAAADAQABAAABAQCetCeLvuaG/Nzc7tuEZ6+eroHWzNg+V7XjzUZZzIiPs35YPFue4guDzEkROi5ThN6/15Nr+lHpoe5yuo4IcCVYyPWSCEc5pymU7r855ilFzfU77RlvBl6sXPVFbg45vEigfpwMTt2e/2VFSAP2lIuli+RatyPmlRVdqBY8lxgiiGSTBFVbDlFTySQyYeyjUI3XcebV1A2NTX637I3LaSQGxSFH3Y8cavFvNYOQS5f7LREzl30+pYd7nh/4jQwQ/UIRfR3tDAUHnlzCUtfcx4ViWRs00Z5nYDr2Z2uzaSjzQvZRcmCMMk/Bp5inHf3q2zNdh4ETIIRqzthpN3eBGGhp RNASeqDevKey" | sudo tee /home/jyin/.ssh/authorized_keys

sudo echo "ssh-rsa AAAAB3NzaC1yc2EAAAADAQABAAABAQCWqH042BJDNeCjWRDbgWbyY4XIroaialZavISW17J8mzOUDyPHMWfK6mkeYOVIIpNYnJ0X1iaLstd7f2Mf7g3WsJChr4gnrQcRY8/kpHQPqbQUK52RfNf33kfEtzIpJiSCxYu9Vjg9ZvhUmJL7dYxdE4ikTjATg6LQBCQGrELqi8KbVQ2HM8jM+C9Q4seUz8p6qaz/EOvrRvwPi22Z0xpNWPzzin2d0jWIWiRIrkMVPGEFWgyqRf3d10F980nAUSkytRiSWYoHVNUgKJOlzCeAlGKSVhAFjSy9CMp8OpXqNKbldHPPgqViRjqckR4yFfanXysyGWLkGSdkepa3W0xZ
ssh-rsa AAAAB3NzaC1yc2EAAAADAQABAAABAQCgtgNraAZilw2X/QZsBSTRBK+fX7Vh5dYCd6Xpps+K6licSHEuHNdFX/5ktxPhvsGvEoopOQYNs+sp738J0XzM/NZ3IzbkL5bGomrK7kF4AgQ6c54u+PgdM39zTD8Jvnxw6C94oTbxM9uHhtfviOmqdr/C1sSTmPV+jyUCQKciWpiNL/to6EG6bRU+vuB23gJ816S21+t6uQXM8+7jZOKbFRNYVueyzX1AN3KaZ9UFnCLl5PLeQLxZKK+3U6afcz74AlmsepBrjDKEZRAvzG5nUh0w+KzqBWgm0dilhS2Sv8V7jVjLdAjNVZackkOzXfREILczeu5VDnphdnZPCxyf" | sudo tee /home/qmeng/.ssh/authorized_keys

sudo echo "ssh-rsa AAAAB3NzaC1yc2EAAAADAQABAAABAQCWqH042BJDNeCjWRDbgWbyY4XIroaialZavISW17J8mzOUDyPHMWfK6mkeYOVIIpNYnJ0X1iaLstd7f2Mf7g3WsJChr4gnrQcRY8/kpHQPqbQUK52RfNf33kfEtzIpJiSCxYu9Vjg9ZvhUmJL7dYxdE4ikTjATg6LQBCQGrELqi8KbVQ2HM8jM+C9Q4seUz8p6qaz/EOvrRvwPi22Z0xpNWPzzin2d0jWIWiRIrkMVPGEFWgyqRf3d10F980nAUSkytRiSWYoHVNUgKJOlzCeAlGKSVhAFjSy9CMp8OpXqNKbldHPPgqViRjqckR4yFfanXysyGWLkGSdkepa3W0xZ" | sudo tee /home/yxu/.ssh/authorized_keys

sudo echo "ssh-rsa AAAAB3NzaC1yc2EAAAADAQABAAABAQCwp9oKUanneTeu4KfKaSVZlWNOlokZ0L6CJ+nGdvxsPfYv0+Jt6PWh65J9Qsr1Nbm6oYJGkOU4V8akA+8q+D2Ebe5sk+6l5uM94+XBVvEzt8KWE7FUNUkdrD93dYWF2HuNqQ5gsfZMN8yTRyDsaF09xyS2PIlMU1EwkWlsvBszfST1kfBrTwwnTcd2f8cybXtzagfH0SuXhYk7c+R/LPBz2qc9pfDa8B3PHIIeoK0UNfoUB0KWXtDi46vC5CA9qHltzOYqRl6UQyFgp9W9XRxG3Vg5FzzMxQ/cYCdtLKLAE4QmTGn1TAHqfyq6oZ6twuP6fIG0MSpA4ufta04HoKIH" | sudo tee /home/cbadger/.ssh/authorized_keys

sudo echo "ssh-rsa AAAAB3NzaC1yc2EAAAADAQABAAABAQCFv7S2F+uqzEUOKmZWIl8TNGzIyWbe2wFaoxn3Ok2g9p5jDXOng0tGKfHdj7muWrF8BsB/hxn/9KfH8kGOsqAoNlDs6kRBBd8RoIrkzOaxiqa8YxIQm907nYfQbGPHpHDGD/N/BePl09FEGVYjxXqBiOI169+JEoYf9s7xKQ6bXte+ZSWQZpcSsBVEtf6AW/OjNqmWigwjUxRYrMo9l3+e/dTEWEemSH4ihYHR4JUJNEIN4RvcxRZdHx5Rco7H6g7EtWttQsMoMi2WBOcqsMw8evDn0znmj8goK+W0gCRiWdRDuGpfHhbWkxblU18xP" | sudo tee /home/jpatterson/.ssh/authorized_keys


for user in jyin qmeng yxu cbadger jpatterson; 
do 
    sudo chown $user.$user /home/$user/.ssh
done