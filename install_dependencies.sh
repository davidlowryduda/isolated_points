#!/bin/bash
currdir=$(pwd)
read -p "Enter the directory to download parallel2023: " whereparallel
read -p "Enter the directory to clone the OpenImage repo: " whereopenimage
 
#Install parallel
    cd $whereparallel
    wget https://ftpmirror.gnu.org/parallel/parallel-20230222.tar.bz2
    wget https://ftpmirror.gnu.org/parallel/parallel-20230222.tar.bz2.sig
    gpg parallel-20230222.tar.bz2.sig
    bzip2 -dc parallel-20230222.tar.bz2 | tar xvf -
    cd parallel-20230222
    ./configure && make && sudo make install
    cd ..

#Clone David Roe's version of the OpenImage directory 
cd $whereopenimage
    git clone https://github.com/roed314/OpenImage
cd $currdir

