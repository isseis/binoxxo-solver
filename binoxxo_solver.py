#! /usr/bin/python
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


'''
Dump the board state.

Args:
    m: Curretn board state.
'''
def dump(m):
    print '  0 1 2 3 4 5 6 7 8 9'
    for i, l in enumerate(m):
        print i, ' '.join(l)
    print ''


'''
Transpose the matrix in place.
'''
def trans(m):
    for i in range(10):
        for j in range(i+1, 10):
            m[i][j], m[j][i] = m[j][i], m[i][j]


def print_debug(m, name, direction, changed):
    if debug and changed:
        print '%s[%s] changed=%d' % (name, direction, changed)
        dump(m)


'''
Apply the given algorithm to each row or column.

Args:
    m: The current board state.
    f: Algorighm to compute the next column (row) state from the current column
       (row) state.
    name: The readable name of the algorithm 'f'.

Returns:
    Number of cells modified.
'''
def step(m, f, name):
    changedH = 0
    changedV = 0

    for l in m:
        changedH += f(l)
    print_debug(m, name, 'H', changedH)
    trans(m)
    for l in m:
        changedV += f(l)
    trans(m)
    print_debug(m, name, 'V', changedV)

    return changedH + changedV


'''
Apply the given algorithm to each row or column pair.

Args:
    m: The current board state.
    f: Algorighm to compute the next column (row) state from the current column
       (row) state.
    name: The readable name of the algorithm 'f'.

Returns:
    Number of cells modified.
'''
def step2(m, f, name):
    changedH = 0
    changedV = 0

    for l1 in m:
        for l2 in m:
            changedH += f(l1, l2)
    print_debug(m, name, 'H', changedH)
    trans(m)
    for l1 in m:
        for l2 in m:
            changedV += f(l1, l2)
    trans(m)
    print_debug(m, name, 'V', changedH)
  
    return changedH + changedV


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
Generates permutation.

Args:
    a: number of 'X' in the result
    b: number of 'O' in the result
'''
def gen(a, b):
    if a == 0 and b == 0:
        return []
    elif a == 0:
        e = []
        for i in range(b):
            e.append('O')
        return [e]
    elif b == 0:
        e = []
        for i in range(a):
            e.append('X')
        return [e]
    else:
        r = []
        for l in gen(a-1, b):
            r.append(['X'] + l)
        for l in gen(a, b-1):
            r.append(['O'] + l)
        return r


'''
Fill characters into empty charcters in the list, and creates a new list.

Args:
    L: The list containing both empty and non-empty characters. The empty characters
       will be replaced by characters in the other list.
    l: The list containing non empty characters to be merged into L. Will be overridden.

Example:
    merge(['O', ' ', ' ', 'X'], ['X', 'O']) -> ['O', 'X', 'O', 'X' ]
'''
def merge(L, l):
    if L.count(' ') != len(l):
        raise RuntimeError('Number of element mismatch. L=' + str(L) + ', l=' + str(l))

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

    for g in gen(5 - l.count('X'), 5 - l.count('O')):
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
def compare(l, L):
    changed = 0
    if l.count(' ') > 0 and L.count(0) == 0:
        changed += compare_sub(l, L, 'X')
        changed += compare_sub(l, L, 'O')
    return changed


def compare_sub(l, L, x):
    Lx = set([i for i, v in enumerate(L) if v == x])
    lx = set([i for i, v in enumerate(l) if v == x])
    if len(Lx) == 5 and len(lx) == 4 and lx.issubset(Lx):
        i = Lx.difference(lx).pop()
        if l[i] == ' ':
            l[i] = flip(x)
            return 1
    return 0


def is_solved(m):
    for l in m:
        if l.count('O') != 5 or l.count('X') != 5:
            return False
    trans(m)
    for l in m:
        if l.count('O') != 5 or l.count('X') != 5:
            return False
    trans(m)
    return True


def validate(m):
    ret = True
    for l in m:
        if not validate_line(l):
            ret = False
    trans(m)
    for l in m:
        if not validate_line(l):
            ret = False
    trans(m)
    return ret


def validate_line(l):
    return (l.count('O') + l.count('X') + l.count(' ') == len(l)
        and validate_line_sub(l, 'O')
        and validate_line_sub(l, 'X'))


def validate_line_sub(l, x):
    if l.count(x) > len(l) / 2:
        return False
    for i in range(len(l) - 2):
        if l[i] == x and l[i+1] == x and l[i+2] == x:
            return False
    return True


def main(b):
    # "OX\n"
    # "XO\n"
    # => [ ['O','X'], ['X','O'] ]
    m = [list(l) for l in b.split('\n')]
    print 'initial'
    dump(m)

    while (step(m, sequence, 'sequence')
            or step(m, middle, 'middle')
            or step2(m, compare, 'compare')
            or step(m, fill_try, 'fill_try')):
        if debug and not validate(m):
            print 'validation failure!'
            dump(m)
            sys.exit(1)

    print 'final'
    dump(m)
    print 'Solved!' if is_solved(m) else 'Failed to solve.'


if __name__ == '__main__':
    main(board)