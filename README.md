# Code Cost Simulations for Genetic Code Error-Minimization Evaluations

This code calculates the cost for the mutation of all codons to each other codon in the codon set. The cumulative code cost of an arrangement of codon-to-amino acid mappings is calcuated. Additionally, a sample size is read in as an argument and is used to generate n (sample size) randoly assigned codon-to-amino acid mappings, and the costs associated with each codon mutation and the cumulative cost of each randomly generated structure of the genetic code is calculated.

Result statistics and plots are generated as output and stored in the specified output folder.

The time taken to run the program is printed to the terminal once the program is finished running.

The amino acid difference values are based off of multiple different models, namely 'Koonin' utilising amino acid polar requirments,
'Higgs' utilising Higgs' multiple amino acid property distance matrix, 'PAM250' utilising the PAM250 matrix, and 'SeqPredNN' utilising the SeqPredNN neural network to determine amino acid similarities.

### Executing program

* Input files: all should be placed in the input-matrices/ folder (see existing files for format specifics)
- the Higgs amino acid distance matrix (.csv)
- SeqPredNN amino acid substitution matrix (.csv)
- PAM250 matrix (.csv)

* Run the following command to execute the code:

```
python3 <path to cost_function.py> <sample size (number of randomly generate codon-to-amino acid mapping simulations)> <mode, ie. standard/primordial> <path to SeqPredNN amino acid substitution matrix> <path to Higgs amino acid distance matrix> <path to PAM250 matrix> <path to output folder>
```

For example:
```
python3 cost_function/cost_function.py 100000 primordial cost_function/input-matrices/oversampled_normalised_matrix.csv cost_function/input-matrices/Higgs_aa_distance_matrix.csv cost_function/input-matrices/pam250.csv cost_function/output
```

## Authors

Lara Berry
23594292@sun.ac.za
