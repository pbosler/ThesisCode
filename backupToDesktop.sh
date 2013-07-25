#!/bin/bash

cd /Volumes/Warehouse/IUTAMPaper/sphereCode

echo 'Entering : ' `pwd`

echo 'Backing up to desktop folder...'

ls *.f90 
ls *.f
ls *.cpp
ls Make*
ls *.namelist
ls *.sh
ls *.ncl

cp *.f90 ~/Desktop/ThesisCode/.
cp *.f ~/Desktop/ThesisCode/.
cp *.cpp ~/Desktop/ThesisCode/.
cp Makefile ~/Desktop/ThesisCode/.
cp *.namelist ~/Desktop/ThesisCode/.
cp *.sh ~/Desktop/ThesisCode/.
cp *.ncl ~/Desktop/ThesisCode/.

TODAY=$(data +"%m_%d_%Y")

tar -czvf ~/Desktop/thesisCodeExport_$TODAY.tar.gz *.f90 *.f Makefile *.cpp *.sh *.namelist