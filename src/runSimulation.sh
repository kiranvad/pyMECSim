#!/bin/bash
chmod u+x ./MECSim

print="Running..: "$1
echo $print

# make a copy of Master.inp into sk file
cp $1 Master.inp # Copy input file as Master

# Run Simulation and plot figures
./MECSim > log.txt




