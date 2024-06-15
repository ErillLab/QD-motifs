import numpy as np
import random
from main import config


def compute_apriori_entropy_GC():
    """
    Compute a priori entropy given GC percentage
    :return: a priori entropy (apr_ent)
    """
    # H_{before}(l) = [-\sum_{S\epsilon\Omega}(f(S)\cdot(log_{2}(f(S)]
    # f(S) is the a priori probability given %GC
    at = 1 - config['GC']
    at /= 2
    gc = config['GC'] / 2
    probabilities = np.array([at, gc, gc, at])  # Estimated probabilities of bases [A, C, G, T]
    # probabilities = [x + config['C'] for x in probabilities]  # Add 1e-10 to eliminate 0s
    apr_ent = -np.sum(probabilities * np.log2(probabilities))
    return apr_ent


apriori_entropy = compute_apriori_entropy_GC()


class Organism:
    # Variables
    __L: int = config['L']  # Length of sites
    __N: int = config['N']  # Number of sites
    __total_IC: int = config['IC']
    __Position_IC: list[float] = np.empty(shape=__L, dtype=float)
    __PROB: list[float] = np.empty(shape=__L, dtype=float)
    __GINI: float = config['GINI']
    __MSE_IC: float = 0
    __Motif: np.ndarray[list[str]] = None  # Sites
    __Compute: bool = config['COMPUTE']  # Flag to control computations

    def __str__(self):
        return f'{self.__Position_IC}'

    # def getters
    def get_l(self):
        return self.__L

    def get_n(self):
        return self.__N

    def get_organism(self):
        return self.__Motif

    def get_total_ic(self):
        return self.__total_IC

    def get_position_ic(self):
        return self.__Position_IC

    def get_ic_flag(self):
        return self.__Compute

    def get_gini(self):
        return self.__GINI

    def get_mse_IC(self):
        return self.__MSE_IC

    # def setters
    def set_motif(self, motif):
        self.__Motif = motif

    def set_probabilities(self, probabilities):
        self.__PROB = probabilities

    def set_total_ic(self, ic):
        self.__total_IC = ic

    def set_position_ic(self, position_ic):  # Sets vector of IC per position and flags IC to not compute again
        self.__Position_IC = position_ic

    def set_GINI(self, gini: float) -> None:
        self.__GINI = gini
        self.set_flag(False)  # To not compute IC & GINI if motif hasn't changed

    def set_flag(self, flag: bool) -> None:  # Might change if more flags are added
        self.__Compute = flag

    def set_mse_IC(self, total_ic: float) -> None:  # Compute and set mean squared error based on Target IC
        """
        Computes organism MSE given target_IC = config['L'] * config['BITS_POSITION']
        :param total_ic: Total IC of the motif
        :return:
        """
        self.__MSE_IC = (total_ic - (config['L'] * config['BITS_POSITION'])) ** 2

    # def class methods
    def generate_motif(self) -> None:
        """
        Generate random binding sites that make up the binding motif of the organism
        :return: motif
        """
        dna_letters = ['A', 'C', 'G', 'T']  # State dna letters
        # Generate n (self.__N) random sequences of length (self.__L) of dna letters
        binding_sites = np.array([[random.choice(dna_letters) for _ in range(self.__L)] for _ in range(self.__N)])
        self.set_motif(binding_sites)

    def generate_max_motifs(self, dna_base) -> None:
        """
               Generate binding sites that make up the binding motif of the organism in the last bin
               :return: motif
               """
        dna_letters = ['A', 'C', 'G', 'T']  # Letras de ADN
        num_letters = len(dna_letters)

        if self.__N % num_letters != 0:
            raise ValueError(
                "N should be multiple of 4")

        # Generate binding sites
        binding_sites = np.empty((self.__N, self.__L), dtype=str)

        for col in range((self.__L - 1)):
            column_letters = dna_letters * (self.__N // num_letters)
            random.shuffle(column_letters)
            for row in range(self.__N):
                binding_sites[row, col] = column_letters[row]

        for row in range(self.__N):
            binding_sites[row, self.__L - 1] = dna_base

        self.set_motif(binding_sites)
    def compute_ic(self) -> None:
        """
        Computes information content of the organism and saves the following variables
        :var: probabilities: list of probabilities of each DNA base in each position of the motif (L)
        :var: total_information_content: total information content of the motif
        :var: position_ic: value of information content in each position of the motif
        :return:
        """
        # I(l) = H_{before}(l)-H_{after}(l)= genomic_bg_entropy-[[-\sum_{S_l\epsilon\Omega}(p(S_l)\cdot(log_{2}(p(S_l)]]
        # p(S_l) is the probability of each base
        total_information_content = 0
        motif_list = self.__Motif.tolist()
        probabilities = []
        position_ic = np.empty(shape=(self.__L,), dtype=float)

        for i in range(self.__L):
            base_counts = np.zeros(4)
            for motif in motif_list:
                # Count apparition of bases per column
                counts = np.array([motif[i].count(base) for base in ['A', 'C', 'G', 'T']])
                base_counts = np.add(base_counts, counts)  # General apparition base count
            base_counts = base_counts + 1   # Add 1 to every position to eliminate 0s (PSEUDO COUNTS)
            # Calculate probability of apparition of every base having in mind we added 1 to every position
            base_probabilities = base_counts / (self.__N + 4)

            probabilities.append(base_probabilities)

            value = (- np.sum(base_probabilities * np.log2(base_probabilities)))
            position_information_content = max(0, apriori_entropy - value)  # Compute IC
            position_ic[i] = position_information_content  # Save IC per position
            total_information_content += position_information_content  # Compute total IC

        # Save values in attributes

        self.set_probabilities(np.array(probabilities))
        self.set_position_ic(position_ic)  # Sets vector of IC per position and flags IC to not compute again
        self.set_total_ic(total_information_content)
        self.set_mse_IC(total_information_content)

    def compute_gini(self) -> None:
        """
        Compute gini given self.__Position_IC and controls any possible error
        :return:
        """
        # Formula used = sum{i=1, n}((2i - (length) - 1)*xi)/(length)*sum{i=1, n}xi

        ics = sorted(self.__Position_IC)  # We sort the IC of every position
        # ics = [x + config['C'] for x in ics]  # Add 1e-10 to eliminate 0s
        index = np.arange(1, self.__L + 1)  # Array of index
        gini = max(0, ((np.sum((2 * index - self.__L - 1) * ics)) / (self.__L * np.sum(ics))))  # Formula applied

        # Gini error control
        assert gini >= 0, ('Gini lower than 0 => ', gini, 'Position ICs =>', self.__Position_IC)
        assert gini <= 1, ('Gini greater than 1 =>', gini, 'Position ICs =>', self.__Position_IC)

        # Save gini value
        self.set_GINI(gini)

    def mutation(self, method: str = 'basic') -> None:
        """
        Applies mutations to the binding motif
        :param method: Mutation method
        :return:
        """
        dna_letters: list[str] = ['A', 'C', 'G', 'T']
        if method == 'basic':   # Change random position
            row = random.randint(0, self.__N - 1)  # Get random row
            col = random.randint(0, self.__L - 1)  # Get random column

            dna_letter = random.choice(dna_letters)  # Get random dna letter

            self.__Motif[row][col] = dna_letter  # Swap previous value of the organism in row, col to the new letter

        elif method == 'additional':  # Swap columns of a motif
            row = random.randint(0, self.__N - 1)  # Get motif
            col1 = random.randint(0, self.__L - 1)  # Get column 1
            col2 = random.randint(0, self.__L - 1)  # Get column 2

            while col1 == col2:  # Control it's not the same column
                col2 = random.randint(0, self.__L - 1)

            self.__Motif[row][col1], self.__Motif[row][col2] = (
                self.__Motif[row][col2], self.__Motif[row][col1])  # Swap columns

    def save_tf_binding_sites(self, position: int) -> str:
        """
        Used to save transcription factor binding sites on a json
        :param position: binding site
        :return:
        """
        binding_site: list[str] = self.__Motif[position].tolist()  # Transform ndarray to list
        string = ''.join(binding_site)  # List to string
        return string
