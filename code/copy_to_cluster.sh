#!/bin/bash

###echo "Hello world"

####cp -r -u code/ /z/edwards_lab1/User/bo1scm/Colombia_rangeLims

####rsync -avr --exclude='./code/stan_files' ./code/ /z/edwards_lab1/User/bo1scm/test

#####--exclude='stan_temp'


#cp -u -r `ls -R | grep "code/" | grep -v "stan_files\|^code/$"` /z/edwards_lab1/User/bo1scm/test/

# copy code
echo "Copying code"
cp -u -r code/ /z/edwards_lab1/User/bo1scm/Colombia_rangeLims

# copy data files
echo "Copying data"
cp -u -r data/ /z/edwards_lab1/User/bo1scm/Colombia_rangeLims

echo "Done"