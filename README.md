# parallel-point-processing


grid.c : Defines the data structures of points and cells. Defines all the cells according to its x, y, z limits. Defines 
a function that takes in a point's id as input and allocates the points to its right cell, and subsequently updating the cell's
container. 
          Things to do : Add the mpi utilities to the same.
          
          
Compilation code : 
mpixlc grid-1.c -o out2.xl -lm

          
Sample run.file : 


#!/bin/sh
#SBATCH --nodes=64
#SBATCH --ntasks=4096
#SBATCH --overcommit

srun -o out64_4096_25600_10k_1.log /gpfs/u/home/GGST/GGSTarun/project/out2.xl 25600 inputfile.txt 10000 Query.txt_o




(
//the arguments for the run file are :
//1: total # input points
//2: inputFile.txt
//3: total # input query points
//4: queryFile.txt
)

