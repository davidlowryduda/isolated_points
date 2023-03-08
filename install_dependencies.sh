#!/bin/bash
currdir=$(pwd)
read -p "Enter the directory to download parallel2023: " whereparallel
read -p "Enter the directory to clone modular curves repo: " wheremodcurves
 
#Install parallel
    cd $whereparallel
    wget https://ftpmirror.gnu.org/parallel/parallel-20230222.tar.bz2
    wget https://ftpmirror.gnu.org/parallel/parallel-20230222.tar.bz2.sig
    gpg parallel-20230222.tar.bz2.sig
    bzip2 -dc parallel-20230222.tar.bz2 | tar xvf -
    cd parallel-20230222
    ./configure && make && sudo make install
    cd ..

#Clone David Roe's github directory on Modular curves
cd $wheremodcurves
    git clone https://github.com/roed314/ModularCurves.git
    cd ModularCurves/equations
    rm -rf Gm-Reduce
    rm -rf OpenImage
    git clone https://github.com/SamSchiavone/Gm-Reduce.git
    git clone https://github.com/AndrewVSutherland/OpenImage.git

cd $currdir

