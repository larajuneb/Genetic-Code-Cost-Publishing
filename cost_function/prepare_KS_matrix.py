import csv
import sys
import pandas as pd


# df = pd.read_csv("input-matrices/Katoh_Suga_matrix.csv")

# df.fillna(value=0, inplace=True)

# df.to_csv("input-matrices/Katoh_Suga_matrix-NO-NULL.csv")

# Katoh_Suga_matrix = list(csv.reader(open("input-matrices/Katoh_Suga_matrix-NO-NULL.csv")))
# Katoh_Suga_matrix[0].pop(1)
# Katoh_Suga_matrix[0].pop(len(Katoh_Suga_matrix[0]) - 1)

Katoh_Suga_matrix_STR = list(csv.reader(open("input-matrices/Katoh_Suga_matrix-CLEAN.csv")))
codons = Katoh_Suga_matrix_STR[0]
Katoh_Suga_matrix_STR.remove(Katoh_Suga_matrix_STR[0]) #remove first row
Katoh_Suga_matrix = [[0 for x in range(64)] for y in range(64)]

for row in range(len(Katoh_Suga_matrix_STR)): #remove first column
    Katoh_Suga_matrix_STR[row].remove(Katoh_Suga_matrix_STR[row][0])

def convert_KS_to_float():
    for i in range(len(Katoh_Suga_matrix_STR)): #row
        for j in range(len(Katoh_Suga_matrix_STR[i])): #col
            Katoh_Suga_matrix[i][j] = float(Katoh_Suga_matrix_STR[i][j])

def normalize_KS_matrix(): #normalise matrix such that each row adds up to 1

    for row in Katoh_Suga_matrix:
        row_sum = sum(row)
        for i in range(len(row)):
            raw_val = row[i]
            if row_sum == 0.0:
                row[i] = 0.0
            else:
                norm_val = raw_val / row_sum
                row[i] = norm_val

convert_KS_to_float()
normalize_KS_matrix()

with open("input-matrices/Katoh_Suga_matrix-NORM.csv", 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerow(codons)
    for i in range(64):
        Katoh_Suga_matrix[i].insert(0, codons[i + 1])
        writer.writerow(Katoh_Suga_matrix[i])



# with open("input-matrices/Katoh_Suga_matrix-CLEAN.csv", 'w', newline='') as file:
#     writer = csv.writer(file)
#     for row in Katoh_Suga_matrix:
#         if row is not Katoh_Suga_matrix[0]:
#             row.pop(0)
#             row.pop(len(row) - 1)
#         writer.writerow(row)
