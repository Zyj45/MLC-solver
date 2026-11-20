import random
import unittest

from mdcrt_adjugate_dual import AdjugateDualMDCRT
from mlc_mdcrt_exact.bidirectional import Z, column, equal_lattices, prepare, verify_solution


class AdjugateDualTests(unittest.TestCase):
    def compare(self, matrices, queries):
        problem = {'kind': 'md', 'matrices': matrices}
        expected = prepare(problem, 'global_tri')
        actual = AdjugateDualMDCRT(matrices)
        self.assertTrue(equal_lattices(expected.B, actual.B))
        for rhs in queries:
            answer = actual.solve(rhs)
            self.assertEqual(answer is None, expected.solve(rhs) is None)
            verify_solution(problem, rhs, answer)
        return actual

    def test_original_example_with_shared_factor(self):
        solver = self.compare([[[2, 0], [2, 6]], [[2, 0], [4, 10]]],
                              [[(1, 5), (1, 1)], [(1, 5), (0, 1)]])
        self.assertEqual(solver.row_count, 2)
        self.assertEqual(solver.factor, 2)

    def test_incompatible_despite_pairwise_intersections(self):
        matrices = [[[2, 0], [0, 1]], [[1, 0], [0, 2]], [[1, 0], [-1, 2]]]
        solver = self.compare(matrices, [[(0, 0), (0, 0), (0, 1)]])
        self.assertIsNone(solver.solve([(0, 0), (0, 0), (0, 1)]))

    def test_no_single_adjugate_row_has_full_order(self):
        solver = self.compare([[[2, 0], [0, -3]]], [[(1, 2)], [(-5, 7)]])
        self.assertEqual(solver.row_count, 2)

    def test_unit_blocks_and_input_validation(self):
        solver = self.compare([[[2, 0], [0, 2]], [[0, 2], [-2, 0]]],
                              [[(1, 2), (3, 4)], [(1, 2), (3, 3)]])
        self.assertEqual(solver.row_count, 0)
        with self.assertRaises(ValueError):
            solver.solve([(1, 2)])
        with self.assertRaises(ValueError):
            AdjugateDualMDCRT([[[1, 2], [2, 4]]])

    def test_random_signed_matrices_and_incompatible_residues(self):
        rng = random.Random(20261111)
        for dimension in (1, 2, 3):
            for trial in range(20):
                matrices = []
                for _ in range(rng.randrange(1, 5)):
                    while True:
                        matrix = [[rng.randint(-4, 4) for _ in range(dimension)] for _ in range(dimension)]
                        if Z(matrix).det():
                            break
                    matrices.append(matrix)
                x = column([rng.randint(-4, 4) for _ in range(dimension)])
                compatible = [list(map(int, (x + Z(matrix) * column(range(dimension))).entries())) for matrix in matrices]
                arbitrary = [[rng.randint(-4, 4) for _ in range(dimension)] for _ in matrices]
                with self.subTest(dimension=dimension, trial=trial):
                    self.compare(matrices, [compatible, arbitrary])


if __name__ == '__main__':
    unittest.main()
