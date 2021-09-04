#! /usr/bin/python

import copy
import unittest
import binoxxo_solver

# helper function to create a board status
def b(*l):
    return [list(s) for s in l]

class TestBinoxxoSolver(unittest.TestCase):

    def assertEqualUnordered(self, expected, actual):
        expected_set = frozenset(tuple([tuple(e) for e in expected]))
        actual_set = frozenset(tuple([tuple(list(e) if isinstance(e, str) else e) for e in actual]))
        self.assertEqual(len(expected), len(expected_set),
                "Duplicated elements in expected: " + str(expected))
        self.assertEqual(len(actual), len(actual_set),
                "Duplicated elements in actual: " + str(actual))
        self.assertEqual(expected_set, actual_set)

    def test_genNull(self):
        self.assertEqual(binoxxo_solver.permutations(0,0), [])

    def test_genOne(self):
        self.assertEqualUnordered(binoxxo_solver.permutations(1,0), ['X'])
        self.assertEqualUnordered(binoxxo_solver.permutations(0,1), ['O'])

    def test_gen20(self):
        self.assertEqualUnordered(binoxxo_solver.permutations(2,0), ['XX'])

    def test_gen02(self):
        self.assertEqualUnordered(binoxxo_solver.permutations(0,2), ['OO'])

    def test_gen11(self):
        self.assertEqualUnordered(binoxxo_solver.permutations(1,1), ['XO', 'OX'])

    def test_gen21(self):
        self.assertEqualUnordered(binoxxo_solver.permutations(2,1), ['XXO', 'XOX', 'OXX'])

    def test_gen22(self):
        self.assertEqualUnordered(binoxxo_solver.permutations(2,2),
                ['XXOO', 'XOXO', 'XOOX', 'OXXO', 'OXOX', 'OOXX'])

    def test_merge1(self):
        self.assertEqual(
                binoxxo_solver.merge([' '], ['O']), ['O'])
        self.assertEqual(
                binoxxo_solver.merge(list(' OO'), list('X')), list('XOO'))
        self.assertEqual(
                binoxxo_solver.merge(list('O O'), list('X')), list('OXO'))

    def test_merge2(self):
        self.assertEqual(
                binoxxo_solver.merge(list('  '), list('OX')), list('OX'))
        self.assertEqual(
                binoxxo_solver.merge(list('  X'), list('OX')), list('OXX'))
        self.assertEqual(
                binoxxo_solver.merge(list(' X '), list('OX')), list('OXX'))
        self.assertEqual(
                binoxxo_solver.merge(list('X  '), list('OX')), list('XOX'))

    def test_merge_unmodified(self):
        # The first list should be unmodified.
        l = list(' XOX')
        binoxxo_solver.merge(l, ['X'])
        self.assertEqual(list(' XOX'), l)

    def test_validaet_dup(self):
        self.assertTrue(binoxxo_solver.validate(b(
            'XOXXOOXOXO',
            'OXOXOXXOOX',
            'XOXOXOOXXO',
            'XOOXOXOOXX',
            'OXXOXOXXOO',
            'OOXXOXOXOX',
            'XXOOXXOOXO',
            'OXOXOOXOXX',
            'XOXOXOXXOO',
            'OXOOXXOXOX')))
        self.assertFalse(binoxxo_solver.validate(b(
            'XOXXOOXOXO',
            'OXOXOXXOOX',
            'XOXOXOOXXO',
            'XOOXOXOOXX',
            'OXXOXOXXOO',
            'OXOXOOXOXX',   # DUP
            'X  OXXO  O',
            'OXOXOOXOXX',   # DUP
            'XOXOXOXXOO',
            'OXOOXXOXOX')))

    def test_validate_notrans(self):
        m = [
            list('XOXXOOXOXO'),
            list('XOOXOXXOOX'),
            list('XOXOXOOXXO'),
            list('XOOXOXOOXX'),
            list('OXXOXOXXOO'),
            list('OOXXOXOXOX'),
            list('XXOOXXOOXO'),
            list('OXOXOOXOXX'),
            list('XOXOXOXXOO'),
            list('OXOOXXOXOX'),
            #     ^four Xs
            ]
        M = copy.deepcopy(m)

        self.assertFalse(binoxxo_solver.validate(m))
        self.assertEqual(M, m)

    def test_validate_line(self):
        self.assertTrue(binoxxo_solver.validate_line(' OO '))

    def test_validate_line_count(self):
        self.assertTrue(binoxxo_solver.validate_line('OXOXXO'))
        self.assertTrue(binoxxo_solver.validate_line('OXXOOX'))

        self.assertFalse(binoxxo_solver.validate_line('OOXOOX'))
        self.assertFalse(binoxxo_solver.validate_line('XXOXXO'))
        self.assertFalse(binoxxo_solver.validate_line('OOOOOO'))

    def test_validate_line_consecutive(self):
        self.assertTrue(binoxxo_solver.validate_line('XOOX'))
        self.assertTrue(binoxxo_solver.validate_line('OXXO'))
        self.assertTrue(binoxxo_solver.validate_line('XOXO'))
        self.assertTrue(binoxxo_solver.validate_line('OXOX'))
        self.assertFalse(binoxxo_solver.validate_line('OOOX'))
        self.assertFalse(binoxxo_solver.validate_line('XOOO'))
        self.assertFalse(binoxxo_solver.validate_line('OOOO'))

    @staticmethod
    def _(l):
        return [list(s) for s in l]

    def test_is_board_derived_empty(self):
        self.assertTrue(binoxxo_solver.is_board_derived([], []))
        self.assertTrue(binoxxo_solver.is_board_derived(
            b('OX', 'XO'),
            b('O ', 'X ')))
        self.assertFalse(binoxxo_solver.is_board_derived(
            b('OX', 'XO'),
            b('O ', 'O ')))


if __name__ == '__main__':
    unittest.main()
