import numpy as np

dna_dict = {'A': 0, 'G': 1, 'C': 2, 'T': 3}

# sequence1 = 'GATC'
# sequence2 = 'CATTG'
# n = 4
# path = 'agct_matrix.csv'
# output_path = 'NW_output.txt'

def is_dna(x):
    for i in x:
        if i not in ['A', 'G', 'C', 'T']:
            return False
    return True

def NeedlemanWunsch(n, submatrix_path, output_filename, sequence1, sequence2, GP = -2):
    if not is_dna(sequence1):
        raise ValueError(f'{sequence1} is not a DNA sequence!')
    if not is_dna(sequence2):
        raise ValueError(f'{sequence2} is not a DNA sequence!')

    substitution_matrix = np.loadtxt(submatrix_path, delimiter=',')
    scoring_matrix = np.zeros((len(sequence2) + 1, len(sequence1) + 1))

    for i in range(len(sequence2) + 1):
        scoring_matrix[i][0] = i * GP

    for j in range(len(sequence1) + 1):
        scoring_matrix[0][j] = j * GP

    UDL = [[[''] for _ in range(len(sequence1) + 1)] for _ in range(len(sequence2) + 1)]

    for i in range(1, len(sequence2) + 1):
        for j in range(1, len(sequence1) + 1):
            H_max = []
            indexes = ['diag', 'up', 'left']
            H_max.append(scoring_matrix[i - 1][j - 1] + substitution_matrix[dna_dict[sequence2[i - 1]]][
                dna_dict[sequence1[j - 1]]])
            H_max.append(scoring_matrix[i - 1][j] + GP)
            H_max.append(scoring_matrix[i][j - 1] + GP)
            scoring_matrix[i][j] = max(H_max)
            UDL[i][j] = [indexes[k] for k in range(len(H_max)) if H_max[k] == max(H_max)]
            UDL[i][0] = 'up'
            UDL[0][j] = 'left'

    alignments = []
    al_seq1 = ''
    al_seq2 = ''
    row = len(sequence2)
    col = len(sequence1)
    def all_alignments(row, col, al_seq1, al_seq2):
        if len(alignments) >= n:
            return
        if row == 0 and col == 0:
            scores = 0
            for k in range(len(al_seq1)):
                if al_seq1[k] == '-' or al_seq2[k] == '-':
                    scores -= 2
                else:
                    r = dna_dict[al_seq2[k]]
                    c = dna_dict[al_seq1[k]]
                    scores += substitution_matrix[r][c]
            alignments.append([al_seq1[::-1], al_seq2[::-1], scores])
            return alignments

        if 'diag' in UDL[row][col] and row > 0 and col > 0:
            all_alignments(row - 1, col - 1, al_seq1 + sequence1[col - 1], al_seq2 + sequence2[row - 1])
        if 'up' in UDL[row][col] and row > 0:
            all_alignments(row - 1, col, al_seq1 + '-', al_seq2 + sequence2[row - 1])
        if 'left' in UDL[row][col] and col > 0:
            all_alignments(row, col - 1, al_seq1 + sequence1[col - 1], al_seq2 + '-')

    all_alignments(row, col, al_seq1, al_seq2)
    f = open(output_filename, 'w')
    for i in range(len(alignments)):
        print(f'Global alignment no. {i + 1}:\n{alignments[i][0]}\n{alignments[i][1]}\nScore: {alignments[i][2]}\n')
        f.writelines(
            f'Global alignment no. {i + 1}:\n{alignments[i][0]}\n{alignments[i][1]}\nScore: {alignments[i][2]}\n\n')
    f.close()

# NeedlemanWunsch(n, path, output_path, sequence1, sequence2, GP = -2)
