# Config.json parameter description

## General parameters
| Parameter | Description                                                                    |
|-----------|--------------------------------------------------------------------------------|
| C         | Small value used to compute certain operations so we don't have divisions by 0 |

## Main parameters
| Parameter                            | Description                                                                        |
|--------------------------------------|------------------------------------------------------------------------------------|
| NUMBER_OF_BINS                       | Number of bins                                                                     |
| NUMBER_OF_INITIAL_POPULATION         | Length of initial organisms generated                                              |
| AVERAGE_FITNESS_ITERATIONS_TO_REVIEW | Number of iterations that we will review the values we have had of average fitness |
| MAX_ITERATIONS                       | Maximum number of iterations until we stop                                         |

## Operation parameters
| Parameter               | Description                                                |
|-------------------------|------------------------------------------------------------|
| BITS_POSITION           | Bits per position                                          |
| MAX_BIN_POPULATION_SIZE | Maximum number of organisms in each bin                    |
| VARIATION               | Maximum value that the average fitness can vary up or down |

## Bin parameters
| Parameter       | Description                                  |
|-----------------|----------------------------------------------|
| RANGE           | Range of values that the bin contains inside |
| POPULATION      | Organisms inside the actual bin              |
| BIN_AVERAGE_IC  | Average information content of the bin       |
| AVERAGE_FITNESS | Average fitness (MSE) of the bin             |

## Organism parameters
| Parameter     | Description                                                         |
|---------------|---------------------------------------------------------------------|
| L             | Length of the organism binding motifs                               |
| N             | Number of binding sites (Depth)                                     |
| IC            | Total IC initial value                                              |
| GINI          | Organism Gini initial value                                         |
| COMPUTE       | Flag to control if organism IC and Gini needs to be computed or not |
| BITS_POSITION | Value of Bits per position used to compute Target IC                |
| GC            | G & C percentage of apparition a priori                             |




