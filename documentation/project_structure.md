# Project structure

## Main
Contains the main loop with the needed function calls to operate, the files that we need to read and the generation of the outputs we want.

## Operations
Acts as a ''middle man''  to call everything we need to operate the algorithm so the main only needs to call the general function

## Bin_class
To control the bin metrics and the population within its range (population inside the bin).

## Organism_Class
The main class that computes every metric we need and represents every organism in the population.