#!/bin/bash

git clone https://github.com/arq5x/bedtools2.git
cd bedtools2
make clean; make all
