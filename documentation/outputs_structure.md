# Outputs Structure

## Name
**The name of the files is composed with:**
* The name of what it is (bins, plot average IC, plot best IC, population, TF binding sites...)
* The Length of the motifs.
* The depth of the population.
* The Bits per Position of the execution.


## Insights of every file

### Bins_X.json
**Contains:** 
* The number of the bin
* Average IC
* Average Fitness
* From each organism inside the bin
  * Total IC
  * Gini
  * Fitness (MSE)
  * Position IC
  * Motifs

### Population_X.json
**Contains:** The list of each organism's motifs 

### TF_binding_sites_X.json
**Contains:** Every transcription factor of every site of the organisms.
Separated by the number of the organism and the number of the Bin.

### Plots
Each plot has in its name what it represents.  
Some of them are images of anomalies found, and they have a name that does not follow the previous structure.