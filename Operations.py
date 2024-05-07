import copy
import numpy as np
import random
import json
import Organism_Class
from main import config

NUMBER_OF_INITIAL_POPULATION = config['NUMBER_OF_INITIAL_POPULATION']
MAX_BIN_POPULATION_SIZE = config['MAX_BIN_POPULATION_SIZE']


def bin_range_assignation(bins):
    """
    Assigns to each bin its range
    :param bins: Containers
    :return:
    """
    partition_size = 1 / len(bins)
    min_value = 0.0
    for i in range(len(bins)):
        max_value = min_value + partition_size
        bins[i].set_bin_range((min_value, max_value))
        min_value = max_value


def generate_population():
    """
    Function that generates the original population
    :return: The list of the initial organisms that make up the global population
    """
    population = []
    for i in range(0, NUMBER_OF_INITIAL_POPULATION):
        organism = Organism_Class.Organism()
        organism.generate_motif()
        population.append(organism)

    population_file = ('./outputs/population_L' + str(config['L']) + '_N' + str(config['N'])
                       + '_BP' + str(config['BITS_POSITION']) + '.json')
    with open(population_file, 'w') as file:
        file.write('Population:\n')
        json.dump([motif.get_organism().tolist() for motif in population], file)
        file.write('\n')
    return population


def compute_information_content(population):
    """
        Calls function compute_ic of every organism that needs to be computed
        :param population: General population
        :return:
        """
    for motif in population:
        if motif.get_ic_flag():
            motif.compute_ic()


def compute_gini(population):
    """
    Calls function compute_gini of every organism that needs to be computed
    :param population: General population
    :return:
    """
    for motif in population:
        if motif.get_ic_flag():
            motif.compute_gini()


def sort_bins(population, bins):
    """
    This function sorts population organisms in its corresponding bin
    :param population: General population
    :param bins: Containers
    :return:
    """
    if len(bins) == 10:     # Sort organism in corresponding bin
        for organism in population:
            value = int(organism.get_gini() * len(bins))
            bins[value].set_population(organism)


def reproduction_process(population):
    """
    This function reproduces the current total population, not having in mind the bin where they come from
    :param population: Global population
    :return: new_population: Population generated from reproduction
    """
    random.shuffle(population)
    new_population = []
    for i in range(0, len(population), 2):
        if i == len(population) - 1:  # Control if odd number
            break
        first_parent = population[i]
        second_parent = population[i + 1]

        first_son = copy.deepcopy(first_parent)  # Copy first parent values
        second_son = copy.deepcopy(second_parent)  # Copy second parent values

        if random.choice([True, False]):  # Random choice to decide if recombination or mutation

            length = ((first_son.get_l()) / 2)
            swap_length = random.randint(1, length)  # Selects swap length M=(1 ... L/2)

            # Select position in org 1 (p1 = 1 ... L-M)
            position_swap_org_1 = random.randint(0, (first_son.get_l() - swap_length) - 1)

            # Select position in org 2 (p2 = 1 ... L-M)
            position_swap_org_2 = random.randint(0, (second_son.get_l() - swap_length) - 1)

            for j in range(0, len(first_son.get_organism())):  # Swap sections in all motifs
                original_first = first_son.get_organism()[j]  # First son motif
                original_second = second_son.get_organism()[j]  # Second son motif

                # Swap sections in an alternative variable to not change values for the second son swap
                new_first = np.concatenate((original_first[:position_swap_org_1],
                                            original_second[position_swap_org_2:(position_swap_org_2 + swap_length)],
                                            original_first[(position_swap_org_1 + swap_length):]))

                new_second = np.concatenate((original_second[:position_swap_org_2],
                                             original_first[position_swap_org_1:(position_swap_org_1 + swap_length)],
                                             original_second[(position_swap_org_2 + swap_length):]))

                # Reassign values to the original ones
                first_son.get_organism()[j] = new_first
                second_son.get_organism()[j] = new_second

        else:

            # Apply mutation to parents deep copy
            if random.choice([True, False]):
                # Apply basic mutation
                first_son.mutation()
                second_son.mutation()
            else:
                # Apply additional mutation
                first_son.mutation('additional')
                second_son.mutation('additional')

        # Set flag to true to compute IC again
        first_son.set_flag(True)
        second_son.set_flag(True)

        # Add sons to the new population generated
        new_population.append(first_son)
        new_population.append(second_son)

    return new_population


def control_overpopulation(population, bins):
    """
    Checks if the bin exceeds maximum population and if so it calls selection_process of the bin
    :param population: Global population
    :param bins: Containers
    :return:
    """
    for bin in bins:
        if bin.get_population_size() > MAX_BIN_POPULATION_SIZE:  # Control if overpopulation happens
            # Compute number of competitions
            number_of_competitions = bin.get_population_size() - MAX_BIN_POPULATION_SIZE
            bin.selection_process(population, number_of_competitions)  # Select organisms to stay


def check_stop(bins, bins_averages):
    """
    Checks if we need to stop the execution (Stopping criterion)
    :param bins: Containers
    :param bins_averages: Container Average fitness of the last config['AVERAGE_FITNESS_ITERATIONS_TO_REVIEW'] iterations
    :return: Boolean to control if we stop or not
    """
    check = np.full((len(bins)), False)
    for i, bin in enumerate(bins):
        if len(bin.get_population()) == MAX_BIN_POPULATION_SIZE:
            bin.average_IC()
            bins_averages[i] = np.append(bins_averages[i],
                                         bin.get_average_ic())[-(config['AVERAGE_FITNESS_ITERATIONS_TO_REVIEW']):]
            greater_than_min = all(n >= (bins_averages[i][-1] - config['VARIATION']) for n in bins_averages[i])
            less_than_max = all(n <= (bins_averages[i][-1] + config['VARIATION']) for n in bins_averages[i])
            if greater_than_min and less_than_max:
                check[i] = True
        else:
            return False
    if np.all(check):
        print('ENDING')
    return np.all(check)
