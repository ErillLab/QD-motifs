import random
import numpy as np
from main import config


class Bin:
    # Variables
    __Bin_Range = config['RANGE']  # Interval that represents [0..0,25] [0,25..0,50]...
    __Population = config['POPULATION'].copy()  # Motifs inside this bin
    __Average_Fitness = config['AVERAGE_FITNESS']
    __Average_IC = config['BIN_AVERAGE_IC']

    # def constructor
    def __init__(self, bin_range=__Bin_Range):
        self.__Bin_Range = bin_range
        self.__Population = config['POPULATION'].copy()

    def __str__(self):
        print('Number organisms: ' + str(len(self.__Population)) + ', Average IC: ' + str(self.__Average_IC))
        return f'{len(self.__Population)}, {self.__Average_IC}'     # Return number of organisms and average IC

    # def getters
    def get_bin_range(self):
        return self.__Bin_Range

    def get_population_size(self):
        return len(self.__Population)

    def get_population(self):
        return self.__Population

    def get_average_fitness(self):
        return self.__Average_Fitness

    def get_average_ic(self):
        return self.__Average_IC

    # def setters
    def set_bin_range(self, bin_range):
        self.__Bin_Range = bin_range

    def set_max_population_size(self, population_size):  # If we want to edit max population of the bin
        self.__Population = population_size

    def set_population(self, organism):
        self.__Population.append(organism)

    # def class methods
    def average_fitness(self):
        """
        Compute average fitness of the bin
        :return:
        """
        if len(self.__Population) == 0:
            self.__Average_Fitness = 0.0
        else:
            average_fitness = 0
            for organism in self.__Population:
                average_fitness += organism.get_mse_IC()
            average_fitness /= len(self.__Population)   # Compute overall average fitness
            self.__Average_Fitness = average_fitness    # Save average fitness

    def average_IC(self):
        """
        Compute average IC of the bin
        :return:
        """
        average_ic = 0
        if len(self.__Population) == 0:
            self.__Average_IC = 0.0
        else:
            for organism in self.__Population:
                average_ic += organism.get_total_ic()
            average_ic /= len(self.__Population)
            self.__Average_IC = average_ic

    def selection_process(self, population, number_of_competitions):  # Compete some organisms
        """
        Function to control overpopulation, it removes a random organism if it has an MSE greater
            than another random organism in the bin
        :param population: General population
        :param number_of_competitions: number of competitions that we need to do in order to control overpopulation
        :return:
        """
        for i in range(number_of_competitions):
            first_organism = random.choice(self.__Population)  # Get random organism
            second_organism = random.choice(self.__Population)  # Get random organism

            while first_organism == second_organism:  # Control second organism is not the same as first one
                second_organism = random.choice(self.__Population)

            # If first organism has lower MSE second organism is removed
            if first_organism.get_mse_IC() <= second_organism.get_mse_IC():
                self.__Population.remove(second_organism)  # Delete organism from bin population
                population.remove(second_organism)  # Delete organism from general population
            else:
                self.__Population.remove(first_organism)  # Delete organism from bin population
                population.remove(first_organism)  # Delete organism from general population
