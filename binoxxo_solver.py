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

import copy

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


'''
Dump the board state.

Args:
    m: the curretn board state.
    M: the previous board state.

    In case the previous board state is given, the diff from the previous
    state to the current state is printed as '*' (empty to O) and '+'
    (empty to X).
'''
def dump(m, M=None):
    if M == None:
        M = m
    print '  0 1 2 3 4 5 6 7 8 9'
    for i, (l, L) in enumerate(zip(m, M)):
        print i,
        for v, V in zip(l, L):
            if v == V:
                print v,
            else:
                print '*' if v == 'O' else '+',
        print ''
    print ''


def trans(m):
    for i in range(10):
        for j in range(i+1, 10):
            m[i][j], m[j][i] = m[j][i], m[i][j]


'''
Apply the given function to each row or column.
'''
def step(m, f, name):
    M = copy.deepcopy(m)
    changed = 0
    direction = ''

    for l in m:
        changed += f(l)
        direction = 'H'
    if changed == 0:
        trans(m)
        for l in m:
            changed += f(l)
            direction = 'V'
        trans(m)

    if changed:
        print '%s[%s] changed=%d' % (name, direction, changed)
        dump(m, M)
    return changed


'''
Apply the given function to each row pair or column pair.
'''
def step2(m, f, name):
    M = copy.deepcopy(m)
    changed = 0
    direction = ''

    for l in m:
        for L in m:
            changed += f(l, L)
            direction = 'H'
    if changed == 0:
        trans(m)
        for l in m:
            for L in m:
                changed += f(l, L)
                direction = 'V'
        trans(m)

    if changed:
        print '%s[%s] changed=%d' % (name, direction, changed)
        dump(m, M)
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
e.g. O_O -> OXO
'''
def middle(l):
    L = copy.deepcopy(l)
    changed = 0
    for i in range(8):
        if L[i+1] == ' ' and L[i] != ' ' and L[i] == L[i+2]:
            l[i+1] = flip(L[i])
            changed += 1
    return changed


'''
Put X next to OO
e.g. XOO_ -> XOOX
'''
def sequence(l):
    L = copy.deepcopy(l)
    changed = 0
    for i in range(8):
        if L[i] == ' ' and L[i+1] != ' ' and L[i+1] == l[i+2]:
            changed += 1
            l[i] = flip(L[i+1])
        if L[i] != ' ' and L[i] == L[i+1] and L[i+2] == ' ':
            changed += 1
            l[i+2] = flip(L[i])
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
In case four X are alredy filled, fill O where X must not be filled.

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


def fill_except(l, v, except1, except2):
    changed = 0
    for i in range(10):
        if l[i] == ' ' and i != except1 and i != except2:
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
    changed = 0
    y = flip(x)

    lx = set([i for i, v in enumerate(l) if v == x])
    Lx = set([i for i, v in enumerate(L) if v == x])
    if len(lx) == 4 and lx.issubset(Lx):
        for i in Lx.difference(lx):
            if l[i] == ' ':
                changed += 1
                l[i] = y

    return changed


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


def main(b):
    # "OX\n"
    # "XO\n"
    # => [ ['O','X'], ['X','O'] ]
    m = [list(l) for l in b.split('\n')]
    print 'initial'
    dump(m)

    while (step(m, sequence, 'sequence')
            or step(m, fill5, 'fill5')
            or step(m, fill4, 'fill4')
            or step(m, middle, 'middle')
            or step2(m, compare, 'compare')):
        pass

    print 'final'
    dump(m)
    print 'Solved!' if is_solved(m) else 'Failed to solve.'


main(board)
