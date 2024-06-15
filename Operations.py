from __future__ import annotations

import copy
import numpy as np
import random
import json

from numpy import bool_

from Bin_Class import Bin
from Organism_Class import Organism
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


def generate_population(flag_max_bin: bool = False) -> list[Organism]:
    """
    Function that generates the original population
    :param flag_max_bin: flag to control if population in the bin with range 0.9-1 is given
    :return: The list of the initial organisms that make up the global population
    """
    population = []
    for i in range(0, NUMBER_OF_INITIAL_POPULATION):
        organism = Organism()
        organism.generate_motif()
        population.append(organism)

    if flag_max_bin:
        dna_letters = ['A', 'C', 'G', 'T']
        for base in dna_letters:
            organism = Organism()
            organism.generate_max_motifs(base)
            population.append(organism)


    # population_file = ('./outputs/population_L' + str(config['L']) + '_N' + str(config['N'])
    #                    + '_BP' + str(config['BITS_POSITION']) + '.json')
    # with open(population_file, 'w') as file:
    #     file.write('Population:\n')
    #     json.dump([motif.get_organism().tolist() for motif in population], file)
    #     file.write('\n')
    return population


def compute_information_content(population: list[Organism]):
    """
        Calls function compute_ic of every organism that needs to be computed
        :param population: General population
        :return:
        """
    for motif in population:
        if motif.get_ic_flag():
            motif.compute_ic()


def compute_gini(population: list[Organism]) -> None:
    """
    Calls function compute_gini of every organism that needs to be computed
    :param population: General population
    :return:
    """
    for motif in population:
        if motif.get_ic_flag():
            motif.compute_gini()


def sort_bins(population: list[Organism], bins: list[Bin]) -> None:
    """
    This function sorts population organisms in its corresponding bin
    :param population: General population
    :param bins: Containers
    :return:
    """
    for organism in population:
        value = int(organism.get_gini() * len(bins))
        if value == config['NUMBER_OF_BINS']:
            value = value - 1
        bins[value].set_population(organism)


def reproduction_process(population: list[Organism]) -> list[Organism]:
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
            swap_length = random.randint(1, int(length))  # Selects swap length M=(1 ... L/2)

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


def control_overpopulation(population: list[Organism], bins: list[Bin]) -> None:
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


def check_stop(bins: list[Bin], bins_minimum_fitness: np.array(list[float])) -> bool | bool_ | bool_:
    """
    Checks if we need to stop the execution (Stopping criterion)
    :param bins_minimum_fitness: Container of Minimum Fitness in bin for the last config A.F.I.T.R. iterations
    :param bins: Containers
    :return: Boolean to control if we stop or not
    """
    check = np.full((len(bins)), False)
    for i, bin in enumerate(bins):
        bin.compute_minimum_fitness()
        bins_minimum_fitness[i] = np.append(bins_minimum_fitness[i],
                                            bin.get_minimum_fitness())[
                                  -(config['AVERAGE_FITNESS_ITERATIONS_TO_REVIEW']):]
        # Check no iteration has changed more than the variation of the last value computed downwards
        greater_than_min = all(
            n >= (bins_minimum_fitness[i][-1] - config['VARIATION']) for n in bins_minimum_fitness[i])
        # Check no iteration has changed more than the variation of the last value computed upwards
        less_than_max = all(n <= (bins_minimum_fitness[i][-1] + config['VARIATION']) for n in bins_minimum_fitness[i])
        none_zero = all(n != 0.0 for n in bins_minimum_fitness[i])  # Check no bin maximum is = 0
        if greater_than_min and less_than_max and none_zero:
            check[i] = True
    if np.all(check):
        print('ENDING')
        return True
    return False
