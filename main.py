import json
with open('./config_values/config.json') as file:
    config = json.load(file)

import Operations
import Bin_Class
import matplotlib.pyplot as plt
import numpy as np


NUMBER_OF_BINS = config['NUMBER_OF_BINS']


def main():

    # Bin generation
    bins = [Bin_Class.Bin() for _ in range(NUMBER_OF_BINS)]     # Population

    # Matrix to save the last AVERAGE_FITNESS_ITERATIONS_TO_REVIEW values:
    bins_averages = np.zeros((NUMBER_OF_BINS, config['AVERAGE_FITNESS_ITERATIONS_TO_REVIEW']))

    Operations.bin_range_assignation(bins)

    # First population generation
    population = Operations.generate_population()
    Operations.compute_information_content(population)  # Compute IC
    Operations.compute_gini(population)  # Compute Gini coefficient
    Operations.sort_bins(population, bins)  # Sort ONLY first population in its corresponding bin

    # Reproduction
    i = 0
    while (not Operations.check_stop(bins, bins_averages)) and i <= config['MAX_ITERATIONS']:
        new_population = Operations.reproduction_process(population)    # Reproduce existing population

        Operations.compute_information_content(new_population)  # Compute IC
        Operations.compute_gini(new_population)     # Compute Gini coefficient
        Operations.sort_bins(new_population, bins)  # Sort ONLY new population in its corresponding bin

        population = population + new_population  # Add new population to existing population

        Operations.control_overpopulation(population, bins)
        i += 1
        print(i)

    for bin in bins:    # Compute average IC & Fitness
        bin.average_IC()
        bin.average_fitness()

    bins_file = ('outputs/bins_L' + str(config['L']) + '_N' + str(config['N'])
                 + '_BP' + str(config['BITS_POSITION']) + '.json')
    with open(bins_file, 'w') as f:
        f.write('Final Bins:\n')
        i = 0
        for k, bin in enumerate(bins):
            f.write('Bin ' + str(k) + ': ')
            f.write('\n')
            f.write('Bin Average IC:' + str(bin.get_average_ic()) + '\n')
            f.write('Bin Average Fitness:' + str(bin.get_average_fitness()) + '\n')
            f.write('\n')
            bin_population = bin.get_population()     # Save all organism values
            for i, organism in enumerate(bin_population, start=1):
                f.write(f'Organism {i} Total IC: {organism.get_total_ic()}\n')
                f.write(f'Organism {i} Position ICs: {organism.get_position_ic()}\n')
                f.write(f'Organism {i} Gini: {organism.get_gini()}\n')
                f.write(f'Organism {i} Motifs: ')
                json.dump([organism.get_organism().tolist()], f)
                f.write('\n')
            f.write('\n')

    # Plot bins # organisms and average IC of each one
    organisms = []
    average_ics = []

    for bin in bins:
        # Returns number of organisms and average IC in a string wit a ',' between
        bin_parameters = bin.__str__().split(',')

        organisms.append(int(bin_parameters[0]))     # First part of the string containing number of organisms in bin
        average_ics.append(float(bin_parameters[1]))   # Second part of the string containing average IC in bin

    # Create plot for number of organisms
    population_plot_file = ('./outputs/plot_number_organisms_L' + str(config['L']) + '_N' + str(config['N'])
                            + '_BP' + str(config['BITS_POSITION']) + '.png')
    plt.bar(range(len(organisms)), organisms, label='# of organisms')
    plt.xlabel('Bins')
    plt.ylabel('Number Organisms')
    plt.title('Number of Organisms within each Bin')
    plt.legend()
    plt.savefig(population_plot_file)
    plt.show()


    # Create plot for average IC
    average_IC_plot_file = ('./outputs/plot_average_IC_L' + str(config['L']) + '_N' + str(config['N'])
                            + '_BP' + str(config['BITS_POSITION']) + '.png')
    plt.bar(range(len(average_ics)), average_ics, label='average IC', alpha=0.5)
    plt.xlabel('Bins')
    plt.ylabel('Average IC')
    plt.title('Average IC within each Bin')
    plt.legend()
    plt.savefig(average_IC_plot_file)
    plt.show()


    # Print bin organisms and information
    organism_tf_binding_sites_files = ('./outputs/tf_binding_sites_L' + str(config['L']) + '_N' + str(config['N'])
                                       + '_BP' + str(config['BITS_POSITION']) + '.json')
    with open(organism_tf_binding_sites_files, 'w') as f:
        for k, bin in enumerate(bins, start=1):
            f.write('Bin ' + str(k) + ':\n')
            for j, organism in enumerate(bin.get_population()):
                organism.__str__()
                f.write('Organism ' + str(j) + ':\n')
                i = 0
                while i < organism.get_n():
                    string = organism.save_tf_binding_sites(i)  # Get organism tf binding sites as string
                    json_string = json.dumps(string)    # Dump string in json
                    f.write(json_string)    # Save json in json file
                    f.write('\n')
                    i += 1
            f.write('\n')

    print('END')


if __name__ == '__main__':
    main()
