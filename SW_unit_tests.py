import unittest
import numpy as np
import os
from SmithWaterman import SmithWaterman

class TestSmithWaterman(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.substitution_matrix = np.array([
            [5, -4, -4, -1],
            [-4, 5, -4, -1],
            [-4, -4, 5, -1],
            [-1, -1, -1, 5]
        ])
        cls.submatrix_path = 'test_agct_matrix1.csv'
        np.savetxt(cls.submatrix_path, cls.substitution_matrix, delimiter=',')

    @classmethod
    def remove_all_test_files(cls):
        os.remove(cls.submatrix_path)

    def test_valid_dna_sequences(self):
        output_path = 'test_output1.txt'
        sequence1 = 'TACGGGCC'
        sequence2 = 'TAGCCCT'
        n = 4

        SmithWaterman(n, self.submatrix_path, output_path, sequence1, sequence2, GP=-2)

        with open(output_path, 'r') as f:
            output = f.read()

        self.assertIn('Local alignment no. 1:', output)
        self.assertIn('TACGGGCC', output)
        self.assertIn('TA---GCC', output)

        os.remove(output_path)

    def test_invalid_dna_sequence(self):
        sequence1 = 'AGCTTGAC'
        sequence2 = 'ATXTACGG'
        n = 4

        with self.assertRaises(ValueError) as context:
            SmithWaterman(n, self.submatrix_path, 'output1.txt', sequence1, sequence2, GP=-2)

        self.assertEqual(str(context.exception), 'ATXTACGG is not a DNA sequence!')

    def test_gap_penalty(self):
        output_path = 'gap_penalty_output1.txt'
        sequence1 = 'CACATAG'
        sequence2 = 'AGATA'
        n = 4

        SmithWaterman(n, self.submatrix_path, output_path, sequence1, sequence2, GP=-1)

        with open(output_path, 'r') as f:
            output = f.read()

        self.assertIn('Local alignment no. 1:', output)
        self.assertIn('AC-ATA', output)
        self.assertIn('A-GATA', output)

        os.remove(output_path)

    def test_multiple_alignments(self):
        output_path = 'test_multiple_alignments_output1.txt'
        sequence1 = 'TACGGGCC'
        sequence2 = 'TAGCCCT'
        n = 4

        SmithWaterman(n, self.submatrix_path, output_path, sequence1, sequence2, GP=-2)

        with open(output_path, 'r') as f:
            output = f.read()

        self.assertIn('Local alignment no. 1:', output)
        self.assertIn('Local alignment no. 2:', output)
        self.assertIn('Local alignment no. 3:', output)

        os.remove(output_path)

    def test_n_stop(self):
        output_path = 'test_n_stop1.txt'
        sequence1 = 'TACGGGCC'
        sequence2 = 'TAGCCCT'
        n = 2

        SmithWaterman(n, self.submatrix_path, output_path, sequence1, sequence2, GP=-2)

        with open(output_path, 'r') as f:
            output = f.read()

        self.assertIn('Local alignment no. 1:', output)
        self.assertIn('Local alignment no. 2:', output)
        self.assertIn('Local alignment no. 3:', output)

        os.remove(output_path)

    def test_score(self):
        output_path = 'test_score1.txt'
        sequence1 = 'GATACAGGCC'
        sequence2 = 'CCCATATG'
        n = 4

        SmithWaterman(n, self.submatrix_path, output_path, sequence1, sequence2, GP=-2)

        with open(output_path, 'r') as f:
            output = f.read()

        self.assertIn('Score: 17.0', output)

        os.remove(output_path)

if __name__ == '__main__':
    unittest.main()
