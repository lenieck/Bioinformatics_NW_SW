import unittest
import numpy as np
import os
from NeedlemanWunsch import NeedlemanWunsch

class TestNeedlemanWunsch(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.substitution_matrix = np.array([
            [5, -4, -4, -1],
            [-4, 5, -4, -1],
            [-4, -4, 5, -1],
            [-1, -1, -1, 5]
        ])
        cls.submatrix_path = 'test_agct_matrix.csv'
        np.savetxt(cls.submatrix_path, cls.substitution_matrix, delimiter=',')

    @classmethod
    def remove_all_test_files(cls):
        os.remove(cls.submatrix_path)

    def test_valid_dna_sequences(self):
        output_path = 'test_output.txt'
        sequence1 = 'AGCTTGAC'
        sequence2 = 'ATATAACGG'
        n = 4

        NeedlemanWunsch(n, self.submatrix_path, output_path, sequence1, sequence2, GP=-2)

        with open(output_path, 'r') as f:
            output = f.read()

        self.assertIn('Global alignment no. 1:', output)
        self.assertIn('A--TAT-AACGG', output)
        self.assertIn('AGCT-TG-AC--', output)

        os.remove(output_path)

    def test_invalid_dna_sequence(self):
        sequence1 = 'AGCTTGAC'
        sequence2 = 'ATXTACGG'
        n = 4

        with self.assertRaises(ValueError) as context:
            NeedlemanWunsch(n, self.submatrix_path, 'output.txt', sequence1, sequence2, GP=-2)

        self.assertEqual(str(context.exception), 'ATXTACGG is not a DNA sequence!')

    def test_gap_penalty(self):
        output_path = 'gap_penalty_output.txt'
        sequence1 = 'GATTACA'
        sequence2 = 'GCATGC'
        n = 4

        NeedlemanWunsch(n, self.submatrix_path, output_path, sequence1, sequence2, GP=-1)

        with open(output_path, 'r') as f:
            output = f.read()

        self.assertIn('Global alignment no. 1:', output)
        self.assertIn('G-ATTACA', output)
        self.assertIn('GCATG-C-', output)

        os.remove(output_path)

    def test_multiple_alignments(self):
        output_path = 'test_multiple_alignments_output.txt'
        sequence1 = 'GATTACAACT'
        sequence2 = 'GCATGC'
        n = 4

        NeedlemanWunsch(n, self.submatrix_path, output_path, sequence1, sequence2, GP=-2)

        with open(output_path, 'r') as f:
            output = f.read()

        self.assertIn('Global alignment no. 1:', output)
        self.assertIn('Global alignment no. 2:', output)
        self.assertIn('Global alignment no. 3:', output)

        os.remove(output_path)

    def test_n_stop(self):
        output_path = 'test_n_stop.txt'
        sequence1 = 'GATTACAACT'
        sequence2 = 'GCATGC'
        n = 2

        NeedlemanWunsch(n, self.submatrix_path, output_path, sequence1, sequence2, GP=-2)

        with open(output_path, 'r') as f:
            output = f.read()

        self.assertIn('Global alignment no. 1:', output)
        self.assertIn('Global alignment no. 2:', output)
        self.assertIn('Global alignment no. 3:', output)

        os.remove(output_path)

    def test_score(self):
        output_path = 'test_score.txt'
        sequence1 = 'GATC'
        sequence2 = 'CATTG'
        n = 4

        NeedlemanWunsch(n, self.submatrix_path, output_path, sequence1, sequence2, GP=-2)

        with open(output_path, 'r') as f:
            output = f.read()

        self.assertIn('Score: 3.0', output)

        os.remove(output_path)

if __name__ == '__main__':
    unittest.main()
