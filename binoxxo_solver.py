#! /usr/bin/python3
# -*- coding: utf-8 -*-

'''
BINOXXO Solver

Rule:
1. Füllen Sie das Rätselgitter mit den Zeichen O und X vollständig aus.
2. Es dürfen nicht mehr als zwei aufeinanderfolgende X und O in einer
   Zeile oder Spalte vorkommen.
3. Pro Zeile und Spalte hat es gleich viele X und O.
4. Alle Zeilen und alle Spalten sind einzigartig.

https://www.marktindex.ch/raetsel/binoxxo/
'''

import argparse
import copy
import itertools
import sys

# Initial board status
# 10x10, each cell must be eitehr ' ', 'X' or 'O'.
board = '''
     XX X 
      X XO
 O        
 O  XX   X
       X X
 X   X    
        OX
       X  
   X      
    O     
'''[1:-1]

# Show debug output and validation
debug = True

class ValidationError(Exception):
    pass


'''
Transpose the matrix in place.
'''
def transpose(m):
    n = len(m)
    for i in range(n):
        for j in range(i+1, n):
            m[i][j], m[j][i] = m[j][i], m[i][j]


def noop(*args):
    pass


def all(f, m, make_iter=iter):
    def _sub(f, m, p):
        try:
            p(m)
            for l in make_iter(m):
                if not f(l):
                    return False
        finally:
            p(m)
        return True

    return _sub(f, m, transpose) and _sub(f, m, noop)


class Algorithm:
    '''
    Create Algorithm instance.

    Args:
        f: Function which receives a line of the matrix, and modify it to
            solve the problem. The function shall return the number of cells
            modified.
        name: Name of the algorithm. Used in debug message.
        make_iter: The the function is applied for each row and column of
            the matrix. You can customize how the iterator is created from
            the matrix by specifying a custom functionh here. This is applied
            for the matrix and the transported matrix.
    '''
    def __init__(self, f, name, make_iter=iter):
        self.f = f
        self.name = name
        self.make_iter = make_iter

    def _apply(self, m, p):
        try:
            p(m)
            changed = 0
            for l in self.make_iter(m):
                changed += self.f(l)
            return changed
        finally:
            p(m)

    '''
    Apply the algorithm to each column in the matrix.

    Args:
        m: The current board state.
        is_trans: Transpose the given matrix when applying the algorithm.

    Returns:
        Number of cells modified.
    '''
    def apply(self, m, is_trans):
        changed = self._apply(m, transpose if is_trans else noop)
        if debug and changed:
            self._print_debug(m, self.name, 'V' if is_trans else 'H', changed)
        return changed

    @staticmethod
    def _print_debug(m, name, direction, changed):
        print('%s[%s] changed=%d' % (name, direction, changed))
        dump(m)
        if not validate(m):
            raise ValidationError


'''
Dump the board state.

Args:
    m: Curretn board state.
'''
def dump(m):
    print('  0 1 2 3 4 5 6 7 8 9')
    for i, l in enumerate(m):
        print(i, ' '.join(l))
    print('')


def flip(n):
    if n == 'O':
        return 'X'
    elif n == 'X':
        return 'O' 
    else:
        return ' '


'''
Put X between two Os

Example:
    O_O -> OXO
'''
def middle(l):
    changed = 0
    for i in range(8):
        if l[i+1] == ' ' and l[i] != ' ' and l[i] == l[i+2]:
            l[i+1] = flip(l[i])
            changed += 1
    return changed


'''
Put X next to OO

Example:
    XOO_ -> XOOX
'''
def sequence(l):
    changed = 0
    for i in range(8):
        if l[i] == ' ' and l[i+1] != ' ' and l[i+1] == l[i+2]:
            changed += 1
            l[i] = flip(l[i+1])
        if l[i] != ' ' and l[i] == l[i+1] and l[i+2] == ' ':
            changed += 1
            l[i+2] = flip(l[i])
    return changed


'''
Generates permutations.

Args:
    a: number of 'X' in the result
    b: number of 'O' in the result

Retrurns:
    Permunations as a list of lists.  For example, permutations(1,2)
    is as below.
    [ ['X', 'O', 'O'], ['O', 'X', 'O'], ['O', 'X', 'X'] ]
'''
def permutations(a, b):
    if a == 0 and b == 0:
        return []
    elif a == 0:
        return [['O'] * b]
    elif b == 0:
        return [['X'] * a]
    else:
        r = []
        for l in permutations(a-1, b):
            r.append(['X'] + l)
        for l in permutations(a, b-1):
            r.append(['O'] + l)
        return r

'''
Fill characters into empty charcters in the list, and creates a new list.

Args:
    L: The list containing both empty and non-empty characters. The empty
       characters will be replaced by characters in the other list.
    l: The list containing non empty characters to be merged into L. Will
       be overridden.
Example:
    merge(['O', ' ', ' ', 'X'], ['X', 'O']) -> ['O', 'X', 'O', 'X' ]
'''
def merge(L, l):
    if L.count(' ') != len(l):
        raise RuntimeError('Number of element mismatch. '
                + 'L=' + str(L) + ', l=' + str(l))
    return [l.pop(0) if E == ' ' else E for E in L]


'''
Generates all possible lines, and see if there is a cell which can only
be occupied by 'X' or 'O'.
'''
def fill_try(l):
    if l.count('O') < 3 and l.count('X') < 3:
        return 0

    # result[i] holds a set of characters which can be placed at the l[i]
    result = []
    for i in range(len(l)):
        result.append(set())

    for g in permutations(5 - l.count('X'), 5 - l.count('O')):
        merged = merge(l, g)
        if validate_line(merged):
            for i, e in enumerate(merged):
                result[i].add(e)

    changed = 0
    for i, e in enumerate(l):
        if len(result[i]) == 1 and l[i] == ' ':
            changed += 1
            l[i] = result[i].pop()

    return changed


'''
In case four X are alreay filled, and there is another fully filled line
L whose four of X posiiton are same as this line, then we cannot X in
the current line where the fifth X in L exist.

Otherwise, the current line and the line L becomes same, and violates
the rule 4.
'''
def compare(i):
    def _sub(l, L, x):
        Lx = set([i for i, v in enumerate(L) if v == x])
        lx = set([i for i, v in enumerate(l) if v == x])
        if len(Lx) == 5 and len(lx) == 4 and lx.issubset(Lx):
            i = Lx.difference(lx).pop()
            if l[i] == ' ':
                l[i] = flip(x)
                return 1
        return 0

    return _sub(i[0], i[1], 'X') + _sub(i[0], i[1], 'O')


def is_solved(m):
    return validate(m) and all(
            lambda l: l.count('O') == l.count('X') and l.count(' ') == 0, m)


'''
Validate the intermediate board status.
'''
def validate(m):
    def _sub(m):
        for l1, l2 in itertools.combinations(m, 2):
            # Rule 4. Alle Zeilen und alle Spalten sind einzigartig.
            if l1.count(' ') <= 1 and l1 == l2:
                return False
        return True

    return (len(m) == 10
            and all(lambda l: len(l) == 10, m)
            and all(lambda l: validate_line(l), m)
            and all(lambda l: l[0].count(' ') > 1 or l[0] != l[1],
                m, lambda m: itertools.combinations(m, 2)))


def validate_line(l):
    def _sub(l):
        for i in range(len(l) - 2):
            # Rule 2. Es dürfen nicht mehr als zwei aufeinanderfolgende
            # X und O in einer Zeile oder Spalte vorkommen.
            if l[i] != ' ' and l[i] == l[i+1] == l[i+2]:
                return False
        return True

    return (l.count('O') + l.count('X') + l.count(' ') == len(l)
            # Rule 3. Pro Zeile und Spalte hat es gleich viele X und O.
            # Note some cells can be empty while solving the puzzle.
            and l.count('O') <= len(l) / 2
            and l.count('X') <= len(l) / 2
            and _sub(l))


'''
Check if any non-empty cell in the original board is modified.

Args:
    m: the current board status.
    M: the original board status.

Returs:
    True if the current board status is derived from the original
    borad status, i.e. all non-empty characters in the original
    board status are kept as is.
'''
def is_board_derived(m, M):
    for l, L in zip(m, M):
        for v, V in zip(l, L):
            if V != ' ' and v != V:
                return False
    return True


def solve(m):
    print('initial')
    dump(m)

    M = copy.deepcopy(m)

    algorithms = [
            Algorithm(sequence, 'sequence'),
            Algorithm(middle, 'middle'),
            Algorithm(compare, 'compare',
                lambda m: itertools.permutations(m, 2)),
            Algorithm(fill_try, 'fill_try'),
            ]
    try:
        while True:
            for a in algorithms:
                if a.apply(m, False) > 0 or a.apply(m, True) > 0:
                    break
            else:
                break
    except ValidationError:
        print('Failed! Invalid board status.')
    else:
        print('final')
        dump(m)

        if not is_board_derived(m, M):
            print('Error. The initial board status was overridden.')
        elif not is_solved(m):
            print('Failed to solve.')
        else:
            print('Solved!')



def main():
    parser = argparse.ArgumentParser(description='Binoxxo puzzle solver.')
    parser.add_argument('--infile', type=argparse.FileType('r'))
    args = parser.parse_args()
    if args.infile:
        b = args.infile.read()
    else:
        b = board

    # "OX\n"
    # "XO\n"
    # => [ ['O','X'], ['X','O'] ]
    m = [list(l) for l in b.split('\n')]
    if not validate(m):
        print('Initial board status is invalid: ' + str(m))
        sys.exit(1)

    solve(m)


if __name__ == '__main__':
    main()
