import random
import csv
from statistics import mean
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import math
import matplotlib
from scipy.stats import spearmanr 
from scipy.stats import ks_2samp
import sys
import os
import time

"""This class calculates the cost for the mutation of all codons to each other codon in the codon set. The cumulative code cost of an arrangement of codon-to-amino
acid mappings is calcuated. Additionally, a sample size is read in as an argument and is used to generate n (sample size) randoly assigned codon-to-amino acid
mappings, and the costs associated with each codon mutation and the cumulative cost of each randomly generated structure of the genetic code is calculated.

Result statistics and plots are generated as output and stored in the specified output folder.

The time taken to run the program is printed to the terminal once the program is finished running.
"""

t0 = time.time()

csv_paths = []
for arg in sys.argv:
    csv_paths.append(arg)
csv_paths.pop(0) #remove "cost_function.py" from list

sample_size = int(csv_paths[0])
mode = csv_paths[1]
if mode == "standard":
    dimension = 64
elif mode == "primordial":
    dimension = 16

#generate SeqPredNN normalised substitution matrix list from .csv 
SeqPredNN_norm_substitution_matrix = list(csv.reader(open(csv_paths[2])))
SeqPredNN_norm_substitution_matrix.remove(SeqPredNN_norm_substitution_matrix[0])

# read Higgs matrix from file ("Higgs_aa_distance_matrix.csv")
Higgs_distance_matrix = list(csv.reader(open(csv_paths[3])))
Higgs_distance_matrix.remove(Higgs_distance_matrix[0])

# read PAM250 matrix from file ("pam250.csv")
PAM250 = list(csv.reader(open(csv_paths[4])))
PAM250.remove(PAM250[0])

Katoh_Suga_matrix = list(csv.reader(open(csv_paths[5])))
Katoh_Suga_matrix.remove(Katoh_Suga_matrix[0])

for row in range(len(Katoh_Suga_matrix)): #remove first column
    Katoh_Suga_matrix[row].remove(Katoh_Suga_matrix[row][0])

#convert Katoh-Suga matrix to float values
for i in range(len(Katoh_Suga_matrix)):
        for j in range(len(Katoh_Suga_matrix[i])):
            Katoh_Suga_matrix[i][j] = float(Katoh_Suga_matrix[i][j])

output_dir = csv_paths[6]

mutation_prob_method = csv_paths[7] # "FH" for Freeland and Hurst weights, "KS" for Katoh-Suga probabilities

if not os.path.exists(output_dir): #check for output folder
    os.mkdir(output_dir)

if not os.path.exists(output_dir + "/matrices"): #check for matrices folder
    os.mkdir(output_dir + "/matrices")

if not os.path.exists(output_dir + "/matrices/" + csv_paths[1]): #check for standard/primordial folder within matrices folder
    os.mkdir(output_dir + "/matrices/" + csv_paths[1])

if not os.path.exists(output_dir + "/plots"): #check for plots folder
    os.mkdir(output_dir + "/plots")

if not os.path.exists(output_dir + "/plots/" + csv_paths[1]): #check for standard/primordial folder within matrices folder
    os.mkdir(output_dir + "/plots/" + csv_paths[1])

if not os.path.exists(output_dir + "/stats"): #check for stats folder
    os.mkdir(output_dir + "/stats")

if not os.path.exists(output_dir + "/stats/" + csv_paths[1]): #check for standard/primordial folder within stats folder
    os.mkdir(output_dir + "/stats/" + csv_paths[1])

for row in range(len(SeqPredNN_norm_substitution_matrix)):
    SeqPredNN_norm_substitution_matrix[row].remove(SeqPredNN_norm_substitution_matrix[row][0])

for row in range(len(Higgs_distance_matrix)):
    Higgs_distance_matrix[row].remove(Higgs_distance_matrix[row][0])
    PAM250[row].remove(PAM250[row][0])


def regenerate_codon_matrix(array):
    """re-populates array with all codons depending on specified mode ("standard" or "primordial")

    Args:
        array (list of str): stores all available codons for this particular mode
    """
    array.clear()
    for first in nucleotides:
        for second in nucleotides:
            if mode == "standard":
                for third in nucleotides:
                    codon = first + second + third
                    array.append(codon)
            elif mode == "primordial":
                codon = first + second
                array.append(codon)

#SeqPredNN codon matrix
SeqPredNN_codon_matrix = [[0 for x in range(dimension)] for y in range(dimension)]
SeqPredNN_codon_matrix_NORM = [[0 for x in range(dimension)] for y in range(dimension)]
SeqPredNN_check = []
SeqPredNN_check_NORM = []
SeqPredNN_sample_code_costs = []
SeqPredNN_sample_code_costs_NORM = []

#Koonin codon matrix
Koonin_codon_matrix = [[0 for x in range(dimension)] for y in range(dimension)]
Koonin_codon_matrix_NORM = [[0 for x in range(dimension)] for y in range(dimension)]
Koonin_check = []
Koonin_check_NORM = []
Koonin_sample_code_costs = []
Koonin_sample_code_costs_NORM = []

#Higgs codon matrix
Higgs_codon_matrix = [[0 for x in range(dimension)] for y in range(dimension)]
Higgs_codon_matrix_NORM = [[0 for x in range(dimension)] for y in range(dimension)]
Higgs_check = []
Higgs_check_NORM = []
Higgs_sample_code_costs = []
Higgs_sample_code_costs_NORM = []

#neutral substitutions values
neutral_codon_matrix = [[0 for x in range(dimension)] for y in range(dimension)]
neutral_codon_matrix_NORM = [[0 for x in range(dimension)] for y in range(dimension)]
neutral_check = []
neutral_check_NORM = []
neutral_sample_code_costs = []
neutral_sample_code_costs_NORM = []

#SeqPredNN amino acid substitutions values only
amino_acid_codon_matrix = [[0 for x in range(dimension)] for y in range(dimension)]
amino_acid_codon_matrix_NORM = [[0 for x in range(dimension)] for y in range(dimension)]
amino_acid_check = []
amino_acid_check_NORM = []
amino_acid_sample_code_costs = []
amino_acid_sample_code_costs_NORM = []

#PAM250 substitutions
PAM250_codon_matrix = [[0 for x in range(dimension)] for y in range(dimension)]
PAM250_codon_matrix_NORM = [[0 for x in range(dimension)] for y in range(dimension)]
PAM250_check = []
PAM250_check_NORM = []
PAM250_sample_code_costs = []
PAM250_sample_code_costs_NORM = []

nucleotides = ['U', 'C', 'A', 'G']
codon_set = []
codon = ""
stop = ['UAA', 'UAG', 'UGA'] #stop codons in standard codon table
disregard = ['UA', 'UG'] #codons to disregard in proposed primordial codon table

#generate all codons in set depending on mode (standard/primordial)
regenerate_codon_matrix(codon_set)

if mode == "standard":
    number_of_codons_per_aa = {"HIS": 2, "ARG": 6, "LYS": 2, "GLN": 2, "GLU": 2, "ASP": 2, "ASN": 2, "GLY": 4, "ALA": 4, "SER": 6, "THR": 4, "PRO": 4, "CYS": 2, "VAL": 4, "ILE": 3, "MET": 1, "LEU": 6, "PHE": 2, "TYR": 2, "TRP": 1, "stop": 3}
    amino_acids = ["HIS", "ARG", "LYS", "GLN", "GLU", "ASP", "ASN", "GLY", "ALA", "SER", "THR", "PRO", "CYS", "VAL", "ILE", "MET", "LEU", "PHE", "TYR", "TRP", "stop"]
    codons_per_aa = {"HIS": ['CAU', 'CAC'], "ARG": ['CGU', 'CGC', 'CGA', 'CGG', 'AGA', 'AGG'], "LYS": ['AAA', 'AAG'], "GLN": ['CAA', 'CAG'], "GLU": ['GAA', 'GAG'], "ASP": ['GAU', 'GAC'], "ASN": ['AAU', 'AAC'], "GLY": ['GGU', 'GGC', 'GGA', 'GGG'], "ALA": ['GCU', 'GCC', 'GCA', 'GCG'], "SER": ['UCU', 'UCC', 'UCA', 'UCG', 'AGU', 'AGC'], "THR": ['ACU', 'ACC', 'ACA', 'ACG'], "PRO": ['CCU', 'CCC', 'CCA', 'CCG'], "CYS": ['UGU', 'UGC'], "VAL": ['GUU', 'GUC', 'GUA', 'GUG'], "ILE": ['AUU', 'AUC', 'AUA'], "MET": ['AUG'], "LEU": ['CUU', 'CUC', 'CUA', 'CUG', 'UUA', 'UUG'], "PHE": ['UUU', 'UUC'], "TYR": ['UAU', 'UAC'], "TRP": ['UGG'], "stop": ['UAA', 'UAG', 'UGA']}
elif mode == "primordial":
    number_of_codons_per_aa = {"HIS": 1, "ARG": 1, "ASP": 1, "ASN": 1, "GLY": 1, "ALA": 1, "SER": 2, "THR": 1, "PRO": 1, "VAL": 1, "ILE": 1, "LEU": 2, "disregard": 2}
    amino_acids = ["HIS", "ARG", "ASP", "ASN", "GLY", "ALA", "SER", "THR", "PRO", "VAL", "ILE", "LEU", "disregard"]
    codons_per_aa = {"HIS": ['CA'], "ARG": ['CG'], "ASP": ['GA'], "ASN": ['AA'], "GLY": ['GG'], "ALA": ['AC'], "SER": ['GC', 'AG'], "THR": ['CC'], "PRO": ['UC'], "VAL": ['GU'], "ILE": ['UU'], "LEU": ['CU', 'AU'], "disregard": ['UA', 'UG']}

Higgs_amino_acid_order = ['PHE','LEU','ILE','MET','VAL','SER','PRO','THR','ALA','TYR','HIS','GLN','ASN','LYS','ASP','GLU','CYS','TRP','ARG','GLY']

# PRIMORDIAL CODON TABLE:
    #_______________________________________
    #       |   U   |   C   |   A   |   G   |
    #-------|-------|-------|-------|-------|
    #   U   |  Ile  |  Pro  |   *   |   *   |
    #-------|-------|-------|-------|-------|
    #   C   |  Leu  |  Thr  |  His  |  Arg  |
    #-------|-------|-------|-------|-------|
    #   A   |  Leu  |  Ala  |  Asn  |  Ser  |
    #-------|-------|-------|-------|-------|
    #   G   |  Val  |  Ser  |  Asp  |  Gly  |
    #_______|_______|_______|_______|_______|

N = 5.7 #using the Freeland Hurst weights below, calculated using calculate_N()
FreelandHurst_mutation_weights = {}
FreelandHurst_mutation_weights["3"] = 1/N
FreelandHurst_mutation_weights["1ts"] = 1/N
FreelandHurst_mutation_weights["1tv"] = 0.5/N
FreelandHurst_mutation_weights["2ts"] = 0.5/N
FreelandHurst_mutation_weights["2tv"] = 0.1/N
FreelandHurst_mutation_weights["otherwise"] = 0 #more than 1 position that differs

polar_requirement_scale = {"HIS": 8.4, "ARG": 9.1, "LYS": 10.1, "GLN": 8.6, "GLU": 12.5, "ASP": 13, "ASN": 10, "GLY": 7.9, "ALA": 7, "SER": 7.5, "THR": 6.6, "PRO": 6.6, "CYS": 11.5, "VAL": 5.6, "ILE": 4.9, "MET": 5.3, "LEU": 4.9, "PHE": 5, "TYR": 5.7, "TRP": 5.3}

all_N = []
all_N_codons = []

purines = ['A', 'G']
pyrimidines = ['C', 'U']

def tv_or_ts(true, mutant):
    """identifies if codon mutation is a transition or transversion

    Args:
        true (str): letter representing original nucleotide ('A', 'T', 'C' or 'G')
        mutant (str): letter representing mutant nucleotide ('A', 'T', 'C' or 'G')

    Returns:
        str: "ts" if mutation is a transition, "tv" if mutation is a transversion
    """
    if (true in purines and mutant in purines) or (true in pyrimidines and mutant in pyrimidines):
        return "ts"
    elif (true in purines and mutant in pyrimidines) or (true in pyrimidines and mutant in purines):
        return "tv"

def calculate_N():
    """calculates N constant for Freeland and Hurst weights
    """
    for true_codon in codon_set:
        sum_for_N = 0
        first = 0
        second = 0
        third = 0
        for mutant_codon in codon_set:
            if true_codon != mutant_codon:
                sum_of_diffs = 0
                first = 0
                second = 0
                third = 0
                if true_codon[0:1] != mutant_codon[0:1]:
                    first = 1
                if true_codon[1:2] != mutant_codon[1:2]:
                    second = 1
                if true_codon[2:3] != mutant_codon[2:3]:
                    third = 1
                sum_of_diffs = first + second + third

                if sum_of_diffs <= 1: #weight for more than one position difference is 0, so don't need to add that to sum_for_N
                    if first == 1: #differ in first position only
                        if tv_or_ts(true_codon[0:1], mutant_codon[0:1]) == "ts":
                            sum_for_N += FreelandHurst_mutation_weights.get("1ts")
                        elif tv_or_ts(true_codon[0:1], mutant_codon[0:1]) == "tv":
                            sum_for_N += FreelandHurst_mutation_weights.get("1tv")
                    elif second == 1: #differ in second position only
                        if tv_or_ts(true_codon[1:2], mutant_codon[1:2]) == "ts":
                            sum_for_N += FreelandHurst_mutation_weights.get("2ts")
                        elif tv_or_ts(true_codon[1:2], mutant_codon[1:2]) == "tv":
                            sum_for_N += FreelandHurst_mutation_weights.get("2tv")
                    elif third == 1: #differ in third position only
                        sum_for_N += FreelandHurst_mutation_weights.get("3")
        all_N.append(sum_for_N)
        all_N_codons.append(true_codon)

#TODO: insert Katoh-Suga matrix calculations
def get_codon_mutation_prob(true_codon, mutant_codon):
    """calculate the probability of a codon mutation using Freeland and Hurst mutation weights

    Args:
        true_codon (str): original codon
        mutant_codon (str): mutant codon

    Returns:
        float or str: if the codons differ by more than one position, '-' is returned; if the codons are the same, '*' is returned; 
        else the Freeland and Hurst codon mutation probability is returned as a float
    """
    if mutation_prob_method == "FH":
        if true_codon != mutant_codon:
            sum_of_diffs = 0
            first = 0
            second = 0
            third = 0
            if true_codon[0:1] != mutant_codon[0:1]:
                first = 1
            if true_codon[1:2] != mutant_codon[1:2]:
                second = 1
            if true_codon[2:3] != mutant_codon[2:3]:
                third = 1
            sum_of_diffs = first + second + third

            if sum_of_diffs <= 1:# differ at only 1 position
                if first == 1: #differ in first position only
                    if tv_or_ts(true_codon[0:1], mutant_codon[0:1]) == "ts":
                        return(FreelandHurst_mutation_weights.get("1ts"))
                    elif tv_or_ts(true_codon[0:1], mutant_codon[0:1]) == "tv":
                        return(FreelandHurst_mutation_weights.get("1tv"))
                elif second == 1: #differ in second position only
                    if tv_or_ts(true_codon[1:2], mutant_codon[1:2]) == "ts":
                        return(FreelandHurst_mutation_weights.get("2ts"))
                    elif tv_or_ts(true_codon[1:2], mutant_codon[1:2]) == "tv":
                        return(FreelandHurst_mutation_weights.get("2tv"))
                elif third == 1: #differ in third position only
                    return(FreelandHurst_mutation_weights.get("3"))
            elif sum_of_diffs > 1: #differ at more than one position
                return '-'
        elif true_codon == mutant_codon:
            return '*'
    elif mutation_prob_method == "KS":
        return Katoh_Suga_matrix[codon_set.index(true_codon)][codon_set.index(mutant_codon)] #THIS LINE PREVENTS AMINO ACID PLOTS

def get_aa_for_codon(codon, codon_dict):
    """accesses amino acid name or stop signal based on specified codon dictionary and codon

    Args:
        codon (str): string specifying codon
        codon_dict (dict of str:list of str): dictionary mapping an amino acid [key] or stop signal to its corresponding assigned codons [value]

    Returns:
        str: name of amino acid or stop signal that codon is assigned to
    """
    for key, value in codon_dict.items():
        if codon in value:
            return key

def get_cost(true_codon_index, mutant_codon_index, model, codon_dict, plot):
    """calculates the cost of a mutation from true codon to mutant codon based on the model being used to quantify the resulting amino acid differences

    Args:
        true_codon_index (int): index of true codon in codon list
        mutant_codon_index (int): index of mutant codon in codon list
        model (str): name of model currently being used for cost function calculations (Higgs, PAM250, SeqPredNN, Koonin, 
        pure SeqPredNN amino acid differences, or neutral amino acid differences)
        codon_dict (dict of str:list of str): dictionary mapping an amino acid [key] or stop signal to its corresponding assigned codons [value]
        plot (bool): states whether or not this calculation is necessary for the production of a plot

    Returns:
        str or float: if the codons differ by more than one position, '-' is returned; if the codons are the same, '*' is returned;
        else the cost of the mutation is returned as a float
    """
    #get codon string name from it's index
    true_codon = codon_set[true_codon_index]
    mutant_codon = codon_set[mutant_codon_index]
    #get amino acid name from codon name
    true_aa = get_aa_for_codon(true_codon, codon_dict)
    mutant_aa = get_aa_for_codon(mutant_codon, codon_dict)

    if true_aa == "stop" or true_aa == "disregard" or mutant_aa == "stop" or mutant_aa == "disregard":
        cost = "-" #cost is disregarded if the mutation if to or from a stop or disregarded codon
    else:
        #get amino acid name from amino acid index
        if model != "Higgs":
            true_aa_index = amino_acids.index(true_aa)
            mutant_aa_index = amino_acids.index(mutant_aa)
        elif model == "Higgs":
            true_aa_index = Higgs_amino_acid_order.index(true_aa)
            mutant_aa_index = Higgs_amino_acid_order.index(mutant_aa)
        elif model == "PAM250":
            true_aa_index = PAM250.index(true_aa)
            mutant_aa_index = PAM250.index(mutant_aa)
        
        # calculate amino acid difference depending on model
        aa_difference = 0
        if model == "SeqPredNN":
            aa_difference = 1 - float(SeqPredNN_norm_substitution_matrix[true_aa_index][mutant_aa_index])
        elif model == "Koonin":
            aa_difference = pow((polar_requirement_scale.get(true_aa) - polar_requirement_scale.get(mutant_aa)), 2)
        elif model == "Higgs":
            aa_difference = Higgs_distance_matrix[true_aa_index][mutant_aa_index]
        elif model == "PAM250":
            aa_difference = float(PAM250[true_aa_index][mutant_aa_index])
            # aa_difference = aa_difference*(-1) #inverted
        elif model == "neutral":
            aa_difference = 1

        # calculate cost of mutation
        codon_mutation_prob = get_codon_mutation_prob(true_codon, mutant_codon)
        if codon_mutation_prob == '*':
            cost = '*'
        elif codon_mutation_prob == '-':
            cost = '-'
        elif model == "amino-acid":
            cost = 1 - float(SeqPredNN_norm_substitution_matrix[true_aa_index][mutant_aa_index]) #THIS COULD BE THE ISSUE
        else:
            cost = float(codon_mutation_prob) * float(aa_difference)
        if plot == True and isinstance(cost, str):
            cost = -10
    
    return (cost)

def normalise_cost(min, max, value):
    """normalises the cost value to be between 0 and 100 based on min and max cost 

    Args:
        min (float): minimum value in the set of costs within which 'value' lies
        max (float): maximum value in the set of costs within which 'value' lies
        value (float): cost value to be normalised between 0 and 100

    Returns:
        float: normalised cost 'value'
    """
    if max - min == 0:
        z = 0
    else:
        z = ((value - min) / (max - min)) * 100
    return(z)

def get_code_cost(cost_array):
    """calculates the overall cost of the code

    Args:
        cost_array (2 dimensional list of floats and str): 2 dimensional list of the costs associated with each possible codon mutation

    Returns:
        float: overall cost of the code
    """
    code_cost = 0
    for row in range(len(cost_array)):
        for cell in range(len(cost_array)):
            if not isinstance(cost_array[row][cell], str):
                code_cost += cost_array[row][cell]
    return code_cost

def generate_sample_set(sample_size, sample_code_costs, sample_code_costs_NORM, model_code_cost, model_code_cost_NORM, model, neutral_cost, neutral_cost_NORM, matrix_length):
    """generates n (sample size) random assignments of codons to amino acids and calculate costs for each structure of the genetic code

    Args:
        sample_size (int): number of times to repeat the generation of random codon-to-amino acid assignments and code cost calculations
        sample_code_costs (list of floats): list of the cumulative code costs for each randomly arranged structure of the genetic code
        sample_code_costs_NORM (list of floats): list of the cumulative code costs for each normalised (0 to 100) randomly arranged structure of the genetic code
        model_code_cost (float): cumulative code cost for the standard genetic code costs calculated using the method specified by the specific model (Higgs, PAM250, 
        SeqPredNN, Koonin, pure SeqPredNN amino acid differences, or neutral amino acid differences)
        model_code_cost_NORM (float): cumulative code cost for the normalised (0 to 100) standard genetic code costs calculated using the method specified by the specific
        model (Higgs, PAM250, SeqPredNN, Koonin, pure SeqPredNN amino acid differences, or neutral amino acid differences)
        model (str): name of model currently being used for cost function calculations (Higgs, PAM250, SeqPredNN, Koonin, 
        pure SeqPredNN amino acid differences, or neutral amino acid differences)
        codon_dict (dict of str:list of str): dictionary mapping an amino acid [key] or stop signal to its corresponding assigned codons [value]
        neutral_cost (float): cumulative code cost for the standard genetic code costs calculated using neutral amino acid differences (all equal to 1)
        neutral_cost_NORM (float): cumulative code cost for the normalised (0 to 100) standard genetic code costs calculated using neutral amino acid differences (all equal to 1)
        matrix_length (int): length of codon matrix, corresponds to number of codons (64 for 'standard' mode, 16 for 'primordial' mode)
    """
    random_codon_assignments = {} #similar to codons_per_aa dict, but instead of true codon assignments, the codons are assigned randomly to amino acids
    leftover = []
    code_costs = []
    norm_code_costs = []

    for i in range(sample_size):
        temp_cost_matrix = [[0 for x in range(matrix_length)] for y in range(matrix_length)]
        norm_cost_matrix = [[0 for x in range(matrix_length)] for y in range(matrix_length)]
        matrix_min_max_check = []
        minimum = 0
        maximum = 0
        regenerate_codon_matrix(leftover)
        for key, value in number_of_codons_per_aa.items():
            random_codon_assignments[key] = [] #a new random codon assignment for each sample, key = amino acid, value = array of codons
            for i in range(number_of_codons_per_aa[key]): #loop through as many times as there are codons for that amino acid
                codon = random.choice(leftover) #choose a random codon for the amino acid/"stop" signal/disregarded codon
                leftover.remove(codon) #remove that codon from the list of available codons
                random_codon_assignments[key].append(codon) #add the randomly chosen codon to the list of codons for that amino acid
        
        for row in range(matrix_length):
            for cell in range(matrix_length):
                temp_cost_matrix[row][cell] = get_cost(row, cell, model, random_codon_assignments, False)

                if not isinstance(temp_cost_matrix[row][cell], str): #add costs for array to check for min and max (excludes "-" for stop/disregarded codons)
                    matrix_min_max_check.append(temp_cost_matrix[row][cell])
        #get code cost for raw data matrix
        minimum = min(matrix_min_max_check)
        maximum = max(matrix_min_max_check)
        code_costs.append(get_code_cost(temp_cost_matrix))

        #get code cost for normalised data matrix
        norm_cost_matrix = normalise_matrix_SAMPLES(minimum, maximum, temp_cost_matrix)

        norm_code_costs.append(get_code_cost(norm_cost_matrix))

        sample_code_costs.append(get_code_cost(temp_cost_matrix))
        sample_code_costs_NORM.append(get_code_cost(norm_cost_matrix))

        random_codon_assignments.clear() #clear random codon assignment for next sample
        temp_cost_matrix = [[0 for x in range(matrix_length)] for y in range(matrix_length)]
        norm_cost_matrix = [[0 for x in range(matrix_length)] for y in range(matrix_length)]
        
    sample_set_stats(code_costs, norm_code_costs, model)

    plot_samples(sample_size, code_costs, model, False, model_code_cost, neutral_cost)
    plot_samples(sample_size, norm_code_costs, model, True, model_code_cost_NORM, neutral_cost_NORM)

    code_costs.clear()
    norm_code_costs.clear()

def sample_set_stats(costs, norm_costs, model):
    """stores statistics summaries for distribution of overall code costs for randomly generated sets codon-to-amino acid assignments

    Args:
        costs (list of floats): list of overall code costs for all randomly generated codon-to-amino acid cost matrices (length = sample size)
        norm_costs (list of floats): list of overall code costs for all randomly generated codon-to-amino acid cost matrices whose cost values have been normalised between 0 and 100
        model (str): name of model currently being used for cost function calculations (Higgs, PAM250, SeqPredNN, Koonin, 
        pure SeqPredNN amino acid differences, or neutral amino acid differences)
    """
    series = pd.Series(costs)
    series_norm = pd.Series(norm_costs)

    filename = output_dir + "/stats/" + mode + "/" + model + "_code_cost_samples_stats.csv"
    with open(filename, mode="w") as file:
        file.write("model, count, mean, std, min, 25%, 50%, 75%, max\n")

        file.write(model + "," + str(series.describe()[0]) + "," + str(series.describe()[1]) + "," + str(series.describe()[2]) + "," + str(series.describe()[3]) + "," + str(series.describe()[4]) + "," + str(series.describe()[5]) + "," + str(series.describe()[6]) + "," + str(series.describe()[7]) + "\n")

        file.write(model + " NORM," + str(series_norm.describe()[0]) + "," + str(series_norm.describe()[1]) + "," + str(series_norm.describe()[2]) + "," + str(series_norm.describe()[3]) + "," + str(series_norm.describe()[4]) + "," + str(series_norm.describe()[5]) + "," + str(series_norm.describe()[6]) + "," + str(series_norm.describe()[7]) + "\n")

def plot_samples(sample_size, costs, model, normalised, code_cost, neutral_cost):
    """plots bar graphs for random samples produced compared to cost of the standard genetic code for that model and for neutral amino acid differences

    Args:
        sample_size (int): number of times to repeat the generation of random codon-to-amino acid assignments and code cost calculations
        costs (list of floats): list of the cumulative code costs for each randomly arranged structure of the genetic code
        model (str): name of model currently being used for cost function calculations (Higgs, PAM250, SeqPredNN, Koonin, 
        pure SeqPredNN amino acid differences, or neutral amino acid differences)
        normalised (bool): specifies if the cost values come from a normalised cost matrix
        code_cost (float): cumulative code cost for the standard genetic code costs calculated using the method specified by the specific model (Higgs, PAM250, 
        SeqPredNN, Koonin, pure SeqPredNN amino acid differences, or neutral amino acid differences)
        neutral_cost (floats): cumulative code cost for the standard genetic code costs calculated using neutral amino acid differences (all equal to 1)
    """
    filename = ""
    title = ""
    costs_temp = []
    for item in costs:
        costs_temp.append(item)
    num_bins = math.isqrt(sample_size)
    mini = float(min(costs))
    maxi = float(max(costs))
    step = (maxi - mini) / num_bins
    bins = {}
    bin_counts = {}
    val = mini
    x = []
    x_bins = []
    y = []
    for i in range(num_bins):
        if i == (num_bins - 1):
            bins[i] = [val, maxi]
        else:
            bins[i] = [val, val + step]
        bin_counts[i] = 0
        val += step

    count = 0
    for key, value in bins.items():
        for item in costs_temp: #not including upper bound, except for in last bin
            if item >= value[0] and item < value[1] and key != (num_bins - 1): 
                bin_counts[key] += 1
            elif item >= value[0] and item <= value[1] and key == (num_bins - 1):
                bin_counts[key] += 1
        x_bins.append("[" + str(round(bins[key][0], 3)) + " - " + str(round(bins[key][1], 3)) + "]")
        x.append(count)
        y.append(bin_counts[key])
        count += 1

    if normalised == False:
        title = "Number of occurrences of randomly simulated code costs per bracket for models generated using the " + model + " method\n Sample size = " + str(sample_size)
        filename = "samples_" + model
    if normalised == True:
        title = "Number of occurrences of randomly simulated code costs per bracket for NORMALISED models generated using the " + model + " method\n Sample size = " + str(sample_size)
        filename = "samples_" + model + "_NORM"    
 
    # create histogram
    middle = (max(y) - min(y))/3

    line_label = model + " standard code cost"
    line_label = line_label.capitalize()

    if model == "amino-acid" and mode == "primordial": #hardcode placement of line label, need to figure this out
    #     diff = max(costs_temp) - min(costs_temp)
    #     lower = min(costs_temp) - 0.1*diff
    #     upper = max(costs_temp) + 0.1*diff
    #     plt.xlim(xmin=lower, xmax=upper)
    #     print("lower: " + str(lower))
    #     print("upper: " + str(upper))
        # print("y-middle: " + str(middle))
        middle = 100
    #     print(y)
    #     print(bins)
    #     print("mini: " + str(mini))
    #     print("maxi: " + str(maxi))
    #     print("step: " + str(step))
    #     print()

    plt.figure(figsize=(15, 10))
    plt.hist(costs, bins=num_bins, edgecolor='white', linewidth=0.3)
    plt.axvline(x = code_cost, color = 'red', linestyle = '--', label = "hi")
    plt.text(code_cost + (0.01*(maxi-mini)), middle, rotation='vertical', s=line_label)

    if model != "neutral" and model != "amino-acid":
        line_label = "Neutral substitutions standard code cost"
        plt.axvline(x = neutral_cost, color = 'green', linestyle = '--', label = line_label)
        plt.text(neutral_cost + (0.01*(maxi-mini)), middle, rotation='vertical', s=line_label)    

    plt.xlabel("Code cost")
    plt.ylabel("Number of occurrences")
    plt.title(title)
    plt.savefig(output_dir + "/plots/" + mode + "/" + f'{filename}.png')
    plt.savefig(output_dir + "/plots/" + mode + "/" + f'{filename}.svg', format="svg")
    plt.close()

def store_cost_matrices(mode, model, matrix):
    """makes csv files of given cost matrix

    Args:
        mode (str): "standard" or "primordial"
        model (str): name of model currently being used for cost function calculations (Higgs, PAM250, SeqPredNN, Koonin, 
        pure SeqPredNN amino acid differences, or neutral amino acid differences)
        matrix (2 dimensional list of floats and str): 2 dimensional list of the costs associated with each possible codon mutation
    """
    filename = output_dir + "/matrices/" + mode + "/" + model + "_cost_matrix.csv"
    codon_string = " ," + ",".join(codon_set) + "\n"
    with open(filename, mode="w") as file:
        file.write(codon_string) #first line
        counter = 0
        for codon_list in matrix:
            file.write(codon_set[counter] + "," + ",".join(map(str, codon_list)) + "\n")
            counter += 1

def calculate_cost_matrix(codon_matrix, check, model):
    """starts matrix cost calculations for this specific model 

    Args:
        codon_matrix (2 dimensional list of floats and str): 2 dimensional list of the costs associated with each possible codon mutation
        check (list of floats): stores all numerical cost values in codon_matrix to later be used to find the minimum and maximum cost values in codon_matrix
        model (str): name of model currently being used for cost function calculations (Higgs, PAM250, SeqPredNN, Koonin, 
        pure SeqPredNN amino acid differences, or neutral amino acid differences)
    """
    for row in range(len(codon_matrix)):
        for cell in range(len(codon_matrix)):
            codon_matrix[row][cell] = get_cost(row, cell, model, codons_per_aa, False)

            if not isinstance(codon_matrix[row][cell], str):
                check.append(codon_matrix[row][cell])

def normalise_matrix(minimum, maximum, original_matrix, norm_check):
    """normalises all values in original_matrix to lie between 0 and 100 

    Args:
        minimum (float): minimum cost value in original matrix
        maximum (float): maximum cost value in original matrix
        original_matrix (2 dimensional list of floats and str): 2 dimensional list of the costs associated with each possible codon mutation
        norm_check (list of floats): stores all numerical cost values in normalised cot matrix to later be used to describe the distribution of these costs

    Returns:
        2 dimensional list of floats and str: storing normalised values from original_matrix (0 to 100)
    """
    normalised = [[0 for x in range(len(original_matrix))] for y in range(len(original_matrix))]

    for row in range(len(original_matrix)):
        for cell in range(len(original_matrix)):
            if not isinstance(original_matrix[row][cell], str):
                normalised[row][cell] = normalise_cost(minimum, maximum, original_matrix[row][cell])
                norm_check.append(normalised[row][cell])
            if isinstance(original_matrix[row][cell], str):
                normalised[row][cell] = original_matrix[row][cell]

    return(normalised)

def normalise_matrix_SAMPLES(minimum, maximum, original_matrix):
    """normalises all values in original_matrix to lie between 0 and 100 

    Args:
        minimum (float): minimum cost value in original matrix
        maximum (float): maximum cost value in original matrix
        original_matrix (2 dimensional list of floats and str): 2 dimensional list of the costs associated with each possible codon mutation

    Returns:
        2 dimensional list of floats and str: storing normalised values from original_matrix (0 to 100)
    """
    normalised = [[0 for x in range(len(original_matrix))] for y in range(len(original_matrix))]

    for row in range(len(original_matrix)):
        for cell in range(len(original_matrix)):
            if not isinstance(original_matrix[row][cell], str):
                normalised[row][cell] = normalise_cost(minimum, maximum, original_matrix[row][cell])
            if isinstance(original_matrix[row][cell], str):
                normalised[row][cell] = original_matrix[row][cell]

    return(normalised)

def spearmans_rank_correlation_tests():
    """performs Spearman's rank correlation tests between all combinations of each model and stores the results in a csv file
    """

    filename = output_dir + "/stats/" + mode + "/Spearmans_rank_correlation_tests.csv"
    with open(filename, mode="w") as file:
        file.write("code cost test, Null hypothesis, Alternative hypothesis, correlation, p-value\n")

        file.write("SeqPredNN vs Koonin, corr == 0, corr != 0," + str(spearmanr(SeqPredNN_check, Koonin_check, alternative='two-sided').correlation) + ","  + str(spearmanr(SeqPredNN_check, Koonin_check, alternative='two-sided').pvalue)+ "\n")
        file.write("SeqPredNN vs Koonin, corr >= 0, corr < 0," + str(spearmanr(SeqPredNN_check, Koonin_check, alternative='less').correlation) + ","  + str(spearmanr(SeqPredNN_check, Koonin_check, alternative='less').pvalue)+ "\n")
        file.write("SeqPredNN vs Koonin, corr <= 0, corr > 0," + str(spearmanr(SeqPredNN_check, Koonin_check, alternative='greater').correlation) + ","  + str(spearmanr(SeqPredNN_check, Koonin_check, alternative='greater').pvalue)+ "\n")

        file.write("SeqPredNN vs Higgs, corr == 0, corr != 0," + str(spearmanr(SeqPredNN_check, Higgs_check, alternative='two-sided').correlation) + ","  + str(spearmanr(SeqPredNN_check, Higgs_check, alternative='two-sided').pvalue)+ "\n")
        file.write("SeqPredNN vs Higgs, corr >= 0, corr < 0," + str(spearmanr(SeqPredNN_check, Higgs_check, alternative='less').correlation) + ","  + str(spearmanr(SeqPredNN_check, Higgs_check, alternative='less').pvalue)+ "\n")
        file.write("SeqPredNN vs Higgs, corr <= 0, corr > 0," + str(spearmanr(SeqPredNN_check, Higgs_check, alternative='greater').correlation) + ","  + str(spearmanr(SeqPredNN_check, Higgs_check, alternative='greater').pvalue)+ "\n")

        file.write("Higgs vs Koonin, corr == 0, corr != 0," + str(spearmanr(Koonin_check, Higgs_check, alternative='two-sided').correlation) + ","  + str(spearmanr(Koonin_check, Higgs_check, alternative='two-sided').pvalue)+ "\n")
        file.write("Higgs vs Koonin, corr >= 0, corr < 0," + str(spearmanr(Koonin_check, Higgs_check, alternative='less').correlation) + ","  + str(spearmanr(Koonin_check, Higgs_check, alternative='less').pvalue)+ "\n")
        file.write("Higgs vs Koonin, corr <= 0, corr > 0," + str(spearmanr(Koonin_check, Higgs_check, alternative='greater').correlation) + ","  + str(spearmanr(Koonin_check, Higgs_check, alternative='greater').pvalue)+ "\n")

def stats():
    """stores a summary of the statistics regarding each model's code costs
    """
    filename = output_dir + "/stats/" + mode + "/code_cost_stats.csv"
    with open(filename, mode="w") as file:
        file.write(" , raw code cost, mean, min, max, NORM code cost, mean NORM, min NORM, max NORM\n")

        #standard
        file.write("SeqPredNN," + str(SeqPredNN_code_cost) + "," + str(mean(SeqPredNN_check)) + "," + str(min(SeqPredNN_check)) + "," + str(max(SeqPredNN_check)) + "," + str(SeqPredNN_code_cost_NORM) + "," + str(mean(SeqPredNN_check_NORM)) + "," + str(min(SeqPredNN_check_NORM)) + "," + str(max(SeqPredNN_check_NORM)) + "\n")

        file.write("Koonin," + str(Koonin_code_cost) + "," + str(mean(Koonin_check)) + "," + str(min(Koonin_check)) + "," + str(max(Koonin_check)) + "," + str(Koonin_code_cost_NORM) + "," + str(mean(Koonin_check_NORM)) + "," + str(min(Koonin_check_NORM)) + "," + str(max(Koonin_check_NORM)) + "\n")

        file.write("Higgs," + str(Higgs_code_cost) + "," + str(mean(Higgs_check)) + "," + str(min(Higgs_check)) + "," + str(max(Higgs_check)) + "," + str(Higgs_code_cost_NORM) + "," + str(mean(Higgs_check_NORM)) + "," + str(min(Higgs_check_NORM)) + "," + str(max(Higgs_check_NORM)) + "\n")

        file.write("Neutral," + str(neutral_code_cost) + "," + str(mean(neutral_check)) + "," + str(min(neutral_check)) + "," + str(max(neutral_check)) + "," + str(neutral_code_cost_NORM) + "," + str(mean(neutral_check_NORM)) + "," + str(min(neutral_check_NORM)) + "," + str(max(neutral_check_NORM)) + "\n")

        file.write("Amino acid," + str(amino_acid_code_cost) + "," + str(mean(amino_acid_check)) + "," + str(min(amino_acid_check)) + "," + str(max(amino_acid_check)) + "," + str(amino_acid_code_cost_NORM) + "," + str(mean(amino_acid_check_NORM)) + "," + str(min(amino_acid_check_NORM)) + "," + str(max(amino_acid_check_NORM)) + "\n")

    filename = output_dir + "/stats/" + mode + "/Kolmogorov_Smirnov_tests.csv"
    with open(filename, mode="w") as file:
        file.write("sample code costs test, statistic, p-value\n")

        file.write("SeqPredNN vs Koonin," + str(ks_2samp(SeqPredNN_sample_code_costs, Koonin_sample_code_costs).statistic) + ","  + str(ks_2samp(SeqPredNN_sample_code_costs, Koonin_sample_code_costs).pvalue)+ "\n")

        file.write("SeqPredNN NORM vs Koonin NORM," + str(ks_2samp(SeqPredNN_sample_code_costs_NORM, Koonin_sample_code_costs_NORM).statistic) + ","  + str(ks_2samp(SeqPredNN_sample_code_costs_NORM, Koonin_sample_code_costs_NORM).pvalue)+ "\n")

        file.write("SeqPredNN vs Higgs," + str(ks_2samp(SeqPredNN_sample_code_costs, Higgs_sample_code_costs).statistic) + ","  + str(ks_2samp(SeqPredNN_sample_code_costs, Higgs_sample_code_costs).pvalue)+ "\n")

        file.write("SeqPredNN NORM vs Higgs NORM," + str(ks_2samp(SeqPredNN_sample_code_costs_NORM, Higgs_sample_code_costs_NORM).statistic) + ","  + str(ks_2samp(SeqPredNN_sample_code_costs_NORM, Higgs_sample_code_costs_NORM).pvalue)+ "\n")

        file.write("Higgs vs Koonin," + str(ks_2samp(Higgs_sample_code_costs, Koonin_sample_code_costs).statistic) + ","  + str(ks_2samp(Higgs_sample_code_costs, Koonin_sample_code_costs).pvalue)+ "\n")

        file.write("Higgs NORM vs Koonin NORM," + str(ks_2samp(Higgs_sample_code_costs_NORM, Koonin_sample_code_costs_NORM).statistic) + ","  + str(ks_2samp(Higgs_sample_code_costs_NORM, Koonin_sample_code_costs_NORM).pvalue)+ "\n")

    SeqPredNN_check.sort()
    Koonin_check.sort()
    Higgs_check.sort()

    #Spearman's rank correlation coefficient:
    spearmans_rank_correlation_tests()

#run calculations for SeqPredNN, Koonin, Higgs, neutral and amino-acid substitution matrices
calculate_cost_matrix(SeqPredNN_codon_matrix, SeqPredNN_check, "SeqPredNN")
calculate_cost_matrix(Koonin_codon_matrix, Koonin_check, "Koonin")
calculate_cost_matrix(Higgs_codon_matrix, Higgs_check, "Higgs")
calculate_cost_matrix(neutral_codon_matrix, neutral_check, "neutral")
calculate_cost_matrix(amino_acid_codon_matrix, amino_acid_check, "amino-acid")
calculate_cost_matrix(PAM250_codon_matrix, PAM250_check, "PAM250")

#create normalised codon mutation cost matrices
SeqPredNN_codon_matrix_NORM = normalise_matrix(min(SeqPredNN_check), max(SeqPredNN_check), SeqPredNN_codon_matrix, SeqPredNN_check_NORM)
Koonin_codon_matrix_NORM = normalise_matrix(min(Koonin_check), max(Koonin_check), Koonin_codon_matrix, Koonin_check_NORM)
Higgs_codon_matrix_NORM = normalise_matrix(min(Higgs_check), max(Higgs_check), Higgs_codon_matrix, Higgs_check_NORM)
neutral_codon_matrix_NORM = normalise_matrix(min(neutral_check), max(neutral_check), neutral_codon_matrix, neutral_check_NORM)
amino_acid_codon_matrix_NORM = normalise_matrix(min(amino_acid_check), max(amino_acid_check), amino_acid_codon_matrix, amino_acid_check_NORM)
PAM250_codon_matrix_NORM = normalise_matrix(min(PAM250_check), max(PAM250_check), PAM250_codon_matrix, PAM250_check_NORM)

#get overall code cost for all matrices
SeqPredNN_code_cost = get_code_cost(SeqPredNN_codon_matrix)
Koonin_code_cost = get_code_cost(Koonin_codon_matrix)
Higgs_code_cost = get_code_cost(Higgs_codon_matrix)
neutral_code_cost = get_code_cost(neutral_codon_matrix)
amino_acid_code_cost = get_code_cost(amino_acid_codon_matrix)
PAM250_code_cost = get_code_cost(PAM250_codon_matrix)

SeqPredNN_code_cost_NORM = get_code_cost(SeqPredNN_codon_matrix_NORM)
Koonin_code_cost_NORM = get_code_cost(Koonin_codon_matrix_NORM)
Higgs_code_cost_NORM = get_code_cost(Higgs_codon_matrix_NORM)
neutral_code_cost_NORM = get_code_cost(neutral_codon_matrix_NORM)
amino_acid_code_cost_NORM = get_code_cost(amino_acid_codon_matrix_NORM)
PAM250_code_cost_NORM = get_code_cost(PAM250_codon_matrix_NORM)

generate_sample_set(sample_size, SeqPredNN_sample_code_costs, SeqPredNN_sample_code_costs_NORM, SeqPredNN_code_cost, SeqPredNN_code_cost_NORM, "SeqPredNN", neutral_code_cost, neutral_code_cost_NORM, dimension)
generate_sample_set(sample_size, Koonin_sample_code_costs, Koonin_sample_code_costs_NORM, Koonin_code_cost, Koonin_code_cost_NORM, "Koonin", neutral_code_cost, neutral_code_cost_NORM, dimension)
generate_sample_set(sample_size, Higgs_sample_code_costs, Higgs_sample_code_costs_NORM, Higgs_code_cost, Higgs_code_cost_NORM, "Higgs", neutral_code_cost, neutral_code_cost_NORM, dimension)
generate_sample_set(sample_size, neutral_sample_code_costs, neutral_sample_code_costs_NORM, neutral_code_cost, neutral_code_cost_NORM, "neutral", neutral_code_cost, neutral_code_cost_NORM, dimension)
generate_sample_set(sample_size, amino_acid_sample_code_costs, amino_acid_sample_code_costs_NORM, amino_acid_code_cost, amino_acid_code_cost_NORM, "amino-acid", neutral_code_cost, neutral_code_cost_NORM, dimension)
generate_sample_set(sample_size, PAM250_sample_code_costs, PAM250_sample_code_costs_NORM, PAM250_code_cost, PAM250_code_cost_NORM, "PAM250", neutral_code_cost, neutral_code_cost_NORM, dimension)

# #generate .csv files
store_cost_matrices(mode, "SeqPredNN", SeqPredNN_codon_matrix)
store_cost_matrices(mode, "SeqPredNN_NORM", SeqPredNN_codon_matrix_NORM)
store_cost_matrices(mode, "Koonin", Koonin_codon_matrix)
store_cost_matrices(mode, "Koonin_NORM", Koonin_codon_matrix_NORM)
store_cost_matrices(mode, "Higgs", Higgs_codon_matrix)
store_cost_matrices(mode, "Higgs_NORM", Higgs_codon_matrix_NORM)
store_cost_matrices(mode, "Neutral_subst", neutral_codon_matrix)
store_cost_matrices(mode, "Neutral_subst_NORM", neutral_codon_matrix_NORM)
store_cost_matrices(mode, "Amino_acid", amino_acid_codon_matrix)
store_cost_matrices(mode, "Amino_acid_NORM", amino_acid_codon_matrix_NORM)
store_cost_matrices(mode, "PAM250", PAM250_codon_matrix)
store_cost_matrices(mode, "PAM250_NORM", PAM250_codon_matrix_NORM)

#generate stats and store in csvs
stats()

t1 = time.time()
total = t1 - t0
print("Total time taken: " + str(total)) 
