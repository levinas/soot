#! /bin/bash

cd programs

git clone https://github.com/samtools/htslib
cd htslib; git pull; make -j; cd ..

git clone https://github.com/samtools/bcftools
cd bcftools; git pull; make -j; cd ..

git clone https://github.com/samtools/tabix
cd tabix; git pull; make -j; cd ..

git clone git://github.com/samtools/samtools.git
cd samtools; git pull; make -j; cd ..

git clone https://github.com/lh3/bwa.git
cd bwa; git pull; make -j; cd ..

git clone https://github.com/arq5x/bedtools
cd bedtools; git pull; make -j; cd ..

git clone --recursive git://github.com/ekg/freebayes.git
cd freebayes; git submodule update --recursive; make -j; cd ..


