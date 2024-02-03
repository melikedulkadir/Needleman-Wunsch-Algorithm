# Melike Nur Dulkadir
# 21992919

import numpy as np
import blosum
import sys


def needleman_wunsch_affine_gap(seq1, seq2, substitution_matrix, gap_open_penalty, gap_extension_penalty):
    # Initialization
    rows = len(seq2) + 1
    cols = len(seq1) + 1

    # Initialize matrices with negative infinity values
    E = np.full((rows, cols), float('-inf'))
    F = np.full((rows, cols), float('-inf'))
    G = np.full((rows, cols), float('-inf'))
    V = np.full((rows, cols), float('-inf'))

    E[0, 0] = float('-inf')
    F[0, 0] = float('-inf')
    G[0, 0] = float('-inf')
    V[0, 0] = 0.0

    # Initialize gap penalties for the first column
    for i in range(1, rows):
        E[i, 0] = gap_open_penalty + (i-1) * gap_extension_penalty
        F[i, 0] = 0.0
        G[i, 0] = 0.0
        V[i, 0] = gap_open_penalty + (i-1) * gap_extension_penalty

    # Initialize gap penalties for the first row
    for j in range(1, cols):
        E[0, j] = 0.0
        F[0, j] = gap_open_penalty + (j-1) * gap_extension_penalty
        G[0, j] = 0.0
        V[0, j] = gap_open_penalty + (j-1) * gap_extension_penalty

    # Fill in the matrices using the Needleman-Wunsch recurrence
    for i in range(1, rows):
        for j in range(1, cols):
            # Calculate E, F, G, and V values. These matrices are initialized and filled with values to represent the
            # optimal alignment scores under different conditions.
            E[i, j] = max(E[i, j - 1] + gap_extension_penalty, V[i, j - 1] + gap_open_penalty)
            F[i, j] = max(F[i - 1, j] + gap_extension_penalty, V[i - 1, j] + gap_open_penalty)
            G[i, j] = V[i - 1, j - 1] + substitution_matrix[(seq1[j - 1])][(seq2[i - 1])]
            V[i, j] = max(G[i, j], E[i, j], F[i, j])

    # Traceback to reconstruct the alignment
    align1, align2 = "", ""
    i, j = rows - 1, cols - 1
    alignment_score = V[i, j]
    while i > 0 and j > 0:
        if V[i, j] == G[i, j]:
            align1 = seq1[j - 1] + align1
            align2 = seq2[i - 1] + align2
            i -= 1
            j -= 1
        elif V[i, j] == E[i, j]:
            align2 = "-" + align2
            align1 = seq1[j - 1] + align1
            j -= 1
        else:
            align2 = seq2[i - 1] + align2
            align1 = "-" + align1
            i -= 1

    return align1, align2, alignment_score


def format_alignment(sequence1, sequence2):     # Format the alignment of two sequences
    alignment = ""
    match_line = ""

    # Create a match line with '|' for matching bases and ' ' for mismatches
    for base1, base2 in zip(sequence1, sequence2):
        if base1 == base2:
            match_line += '|'
        else:
            match_line += ' '

    # Construct the alignment string
    alignment += sequence1 + '\n' + match_line + '\n' + sequence2

    return alignment


def calculate_identity(sequence1, sequence2):       # Calculate sequence identity
    match_count = sum(b1 == b2 for b1, b2 in zip(sequence1, sequence2))
    identity = (match_count / len(sequence1)) * 100
    return str(match_count) + '/' + str(len(sequence1)) + ' (%' + str(identity) + ')'


# Get file path and gap penalties from command-line arguments
file_path = sys.argv[1]
gap_open_penalty = int(sys.argv[2])
gap_extend_penalty = int(sys.argv[3])

with open(file_path, 'r') as input_file:
    sequences = input_file.readlines()

# Read sequences from the input file
sequence1 = sequences[0].strip()
sequence2 = sequences[1].strip()

# The scoring matrix (can be change when the user wants to do so)
# blosum_matrix = blosum.BLOSUM(45)
blosum_matrix = blosum.BLOSUM(62)


# Perform Needleman-Wunsch alignment with affine gap penalties
aligned_seq1, aligned_seq2, alignment_score = needleman_wunsch_affine_gap(sequence1, sequence2, blosum_matrix,
                                                                          gap_open_penalty,
                                                                          gap_extend_penalty)

# Display the alignment results
print(format_alignment(aligned_seq1, aligned_seq2))
print("Alignment score:", alignment_score)
print("Identity value:", calculate_identity(aligned_seq1, aligned_seq2))
