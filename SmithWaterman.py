import numpy as np

dna_dict = {'A': 0, 'G': 1, 'C': 2, 'T': 3}

# sequence1 = 'GATACAGGCC'
# sequence2 = 'CCCATATG'
# n = 100
# path = 'agct_matrix.csv'
# output_path = 'SW_output.txt'

def is_dna(x):
    for i in x:
        if i not in ['A', 'G', 'C', 'T']:
            return False
    return True


def SmithWaterman(n, submatrix_path, output_filename, sequence1, sequence2, GP=-2):
    if not is_dna(sequence1):
        raise ValueError(f'{sequence1} is not a DNA sequence!')
    if not is_dna(sequence2):
        raise ValueError(f'{sequence2} is not a DNA sequence!')

    substitution_matrix = np.loadtxt(submatrix_path, delimiter=',')
    scoring_matrix = np.zeros((len(sequence2) + 1, len(sequence1) + 1))

    max_score = 0
    max_indexes = []

    for i in range(1, len(sequence2) + 1):
        for j in range(1, len(sequence1) + 1):
            diag = scoring_matrix[i - 1][j - 1] + substitution_matrix[dna_dict[sequence2[i - 1]]][
                dna_dict[sequence1[j - 1]]]
            up = scoring_matrix[i - 1][j] + GP
            left = scoring_matrix[i][j - 1] + GP

            scoring_matrix[i][j] = max(0, diag, up, left)

            if scoring_matrix[i][j] > max_score:
                max_score = scoring_matrix[i][j]
                max_indexes = [(i, j)]
            elif scoring_matrix[i][j] == max_score:
                max_indexes.append((i, j))

    alignments = []
    al_seq1 = ''
    al_seq2 = ''

    def all_alignments(row, col, al_seq1, al_seq2):
        if len(alignments) >= n:
            return
        if scoring_matrix[row][col] == 0:
            alignments.append([al_seq1[::-1], al_seq2[::-1], scoring_matrix[max_indexes[0][0]][max_indexes[0][1]]])
            return

        if row > 0 and col > 0 and scoring_matrix[row][col] == scoring_matrix[row - 1][col - 1] + substitution_matrix[
            dna_dict[sequence2[row - 1]]][dna_dict[sequence1[col - 1]]]:
            all_alignments(row - 1, col - 1, al_seq1 + sequence1[col - 1], al_seq2 + sequence2[row - 1])
        if row > 0 and scoring_matrix[row][col] == scoring_matrix[row - 1][col] + GP:
            all_alignments(row - 1, col, al_seq1 + '-', al_seq2 + sequence2[row - 1])
        if col > 0 and scoring_matrix[row][col] == scoring_matrix[row][col - 1] + GP:
            all_alignments(row, col - 1, al_seq1 + sequence1[col - 1], al_seq2 + '-')

    for p in max_indexes:
        all_alignments(p[0], p[1], al_seq1, al_seq2)

    f = open(output_filename, 'w')
    for i in range(len(alignments)):
        print(f'Local alignment no. {i + 1}:\n{alignments[i][0]}\n{alignments[i][1]}\nScore: {alignments[i][2]}\n')
        f.writelines(
            f'Local alignment no. {i + 1}:\n{alignments[i][0]}\n{alignments[i][1]}\nScore: {alignments[i][2]}\n\n')
    f.close()

# SmithWaterman(n, path, output_path, sequence1, sequence2, GP=-2)
