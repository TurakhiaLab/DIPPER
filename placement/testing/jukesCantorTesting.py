import math
import numpy as np

def jukes_cantor_distance(seq1, seq2):
    if len(seq1) != len(seq2):
        raise ValueError("Sequences must be of equal length")
    
    matches = sum(1 for a, b in zip(seq1, seq2) if a == b)
    seq_length = len(seq1)
    p = 1 - (matches / seq_length)
    
    if p >= 0.75:
        return float('inf')
    else:
        return -0.75 * math.log(1 - (4/3) * p)

# Input sequences
sequences = [
    "TTGCGTATGCCATATTCTCATGGATAATTAATGGGGGTTTGC",
    "TTGCGTATGCATTATTCTCATGGATAATTCATGTGGGTTTTC",
    "CTGATTAAGCTAGATTATCATGGATAATTCATATGGGTTGGC",
    "CTGCGTGTGCAATATTCTCATGGATAATTTATGGGGATTTGC"
]

# Number of sequences
n = len(sequences)

# Initialize distance matrix
distance_matrix = np.zeros((n, n))

# Calculate distances
for i in range(n):
    for j in range(i+1, n):
        distance = jukes_cantor_distance(sequences[i], sequences[j])
        distance_matrix[i][j] = distance
        distance_matrix[j][i] = distance  # Matrix is symmetric

# Print the distance matrix
print("Jukes-Cantor Distance Matrix:")
for row in distance_matrix:
    print(" ".join(f"{x:.6f}" for x in row))