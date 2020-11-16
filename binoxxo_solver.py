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
         X
  X    O  
          
O X OO O  
  O       
O    OO X 
 O     X X
  X  X    
  O  X  OO
XO        
'''[1:-1]

# Show debug output and validation
debug = False


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
    changed = 0

    for l in m:
        changed += f(l)
    trans(m)
    for l in m:
        changed += f(l)
    trans(m)

    if debug and changed:
        print '%s changed=%d' % (name, changed)
        dump(m)
    return changed


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
    changed = 0

    for l1 in m:
        for l2 in m:
            changed += f(l1, l2)
    trans(m)
    for l1 in m:
        for l2 in m:
            changed += f(l1, l2)
    trans(m)
  
    if debug and changed:
        print '%s changed=%d' % (name, changed)
        dump(m)
    return changed


def flip(n):
    if n == 'O':
        return 'X'
    elif n == 'X':
        return 'O' 
    else:
        return ' '


'''
Put 'X' between two Os

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
In case five X are already filled, fill O onto the remaining cells.
'''
def fill5(l):
    changed = 0
    changed += fill5sub(l, 'X')
    changed += fill5sub(l, 'O')
    return changed


def fill5sub(l, x):
    changed = 0
    y = flip(x)

    if l.count(x) == 5:
        for i in range(10):
            if l[i] == ' ':
                changed += 1
                l[i] = y

    return changed
 

'''
In case four X are alredy filled, fill O where X can not be filled.

Example:
    __OX_ -> __OXO, O__X_ -> O__XO, _O_X_ -> _O_XO
    Otherwise the fist three cells become OOO, which violates the rule 2
'''
def fill4(l):
    changed = 0
    changed += fill4sub(l, 'X')
    changed += fill4sub(l, 'O')
    return changed


def fill4sub(l, x):
    changed = 0
    y = flip(x)

    if l.count(x) == 4 and l.count(y) <= 3:
        for i in range(8):
            if l[i] == ' ' and l[i+1] == ' ' and l[i+2] == y:
                changed += fill_except(l, y, i, i+1)
            if l[i] == ' ' and l[i+1] == y and l[i+2] == ' ':
                changed += fill_except(l, y, i, i+2)
            if l[i] == y and l[i+1] == ' ' and l[i+2] == ' ':
                changed += fill_except(l, y, i+1, i+2)

    return changed

'''
Fill empty cells by the given value except its index is either e1 or e2.
'''
def fill_except(l, v, e1, e2):
    changed = 0
    for i in range(10):
        if l[i] == ' ' and i != e1 and i != e2:
            changed += 1
            l[i] = v
    return changed


'''
In case four X are alreay filled, and there is another fully filled line
'L' whose four of X posiiton are same as this line, then we cannot X in
the current line where the fifth X in 'L' exist.

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
        l[Lx.difference(lx).pop()] = flip(x)
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
    error = False
    for l in m:
        if l.count('O') > 5 or l.count('X') > 5:
            error = True
    trans(m)
    for l in m:
        if l.count('O') > 5 or l.count('X') > 5:
            error = True
    trans(m)

    if error:
        raise RuntimeError

def main(b):
    # "OX\n"
    # "XO\n"
    # => [ ['O','X'], ['X','O'] ]
    m = [list(l) for l in b.split('\n')]
    print 'initial'
    dump(m)

    try:
        while (step(m, sequence, 'sequence')
                + step(m, fill5, 'fill5')
                + step(m, fill4, 'fill4')
                + step(m, middle, 'middle')
                + step2(m, compare, 'compare')):
            if debug:
                validate(m)
        print 'final'
    except RuntimeError:
        print 'error!'

    dump(m)
    print 'Solved!' if is_solved(m) else 'Failed to solve.'


main(board)
