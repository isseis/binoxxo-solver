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

# Initial board status
#  0: empty
#  1: X
#  2: O
board = [
        [0,0,0,0,0,0,0,0,0,1],
        [0,0,1,0,0,0,0,2,0,0],
        [0,0,0,0,0,0,0,0,0,0],
        [2,0,1,0,2,2,0,2,0,0],
        [0,0,2,0,0,0,0,0,0,0],
        [2,0,0,0,0,2,2,0,1,0],
        [0,2,0,0,0,0,0,1,0,1],
        [0,0,1,0,0,1,0,0,0,0],
        [0,0,2,0,0,1,0,0,2,2],
        [1,2,0,0,0,0,0,0,0,0],
]

def mark(i):
    if i == 1:
        return 'X'
    elif i == 2:
        return 'O'
    else:
        return ' ';


def dump1(l):
    for i in range(len(l)):
        print mark(l[i]),
    print ''

 
def dump():
    print '  0 1 2 3 4 5 6 7 8 9'
    for i in range(len(board)):
        print i,
        dump1(board[i])
    print ''


def main():
    dump()
    while (True):
        if (step(sequence, 'sequence')
                or step(fill5, 'fill5')
                or step(fill4, 'fill4')
                or step(middle, 'middle')
                or step2(compare, 'compare')):
            continue
        break

    solved = True
    for i in range(len(board)):
        if board[i].count(0) > 0:
            solved = False
    if solved:
        print 'Solved!'
    else:
        print 'Failed to solve'


'''
Apply the given function to each line and column.
'''
def step(f, name):
    changed = 0
    for i in range(len(board)):
        changed += f(board[i])

    # create row/colum swapped data, and apply the function
    for c in range(10):
        l = []
        for i in range(10):
            l.append(board[i][c])
        changed += f(l)
        for i in range(10):
            board[i][c] = l[i]

    if changed:
        print '%s changed=%d' % (name, changed)
        dump()
    return changed


'''
Apply the given function to each line pair and column pair.
'''
def step2(f, name):
    changed = 0
    for i in range(len(board)):
        for j in range(len(board)):
            changed += f(board[i], board[j])

    # create row/colum swapped data, and apply the function
    for c in range(10):
        for c2 in range(10):
            l = []
            L = []
            for i in range(10):
                l.append(board[i][c])
                L.append(board[i][c2])
            changed += f(l, L)
            for i in range(10):
                board[i][c] = l[i]

    if changed:
        print '%s changed=%d' % (name, changed)
        dump()
    return changed


def flip(n):
    if n == 1:
        return 2
    elif n == 2:
        return 1
    else:
        return 0

'''
Put 'X' between two Os
e.g. O_O -> OXO
'''
def middle(l):
    changed = 0
    for i in range(len(l) - 2):
        if l[i+1] == 0 and l[i] != 0 and l[i] == l[i+2]:
            l[i+1] = flip(l[i])
            changed += 1
    return changed

'''
Put X next to OO
e.g. XOO_ -> XOOX
'''
def sequence(l):
    changed = 0
    for i in range(len(l)-2):
        if l[i] == 0 and l[i+1] != 0 and l[i+1] == l[i+2]:
            changed += 1
            l[i] = flip(l[i+1])
        if l[i] != 0 and l[i] == l[i+1] and l[i+2] == 0:
            changed += 1
            l[i+2] = flip(l[i])
    return changed


'''
In case five X are already filled, fill O onto the remaining cells.
'''
def fill5(l):
    changed = 0

    if l.count(1) == 5:
        for i in range(len(l)):
            if l[i] == 0:
                changed += 1
                l[i] = 2
    if l.count(2) == 5:
        for i in range(len(l)):
            if l[i] == 0:
                changed += 1
                l[i] = 1
    return changed

def fill_except(l, v, except1, except2):
    changed = 0
    for i in range(len(l)):
        if l[i] == 0 and i != except1 and i != except2:
            changed += 1
            l[i] = v
    return changed

'''
In case four X are alredy filled, fill O where X must not be filled.

Example:
__OX_ -> __OXO, O__X_ -> O__XO, _O_X_ -> _O_XO
Otherwise the fist three cells become OOO, which violates the rule 2
'''
def fill4(l):
    changed = 0

    x = 1
    y = flip(x)

    changed += fill4sub(l, 1)
    changed += fill4sub(l, 2)
    return changed

def fill4sub(l, x):
    changed = 0
    y = flip(x)

    if l.count(x) == 4 and l.count(y) <= 3:
        for i in range(len(l)-2):
            if l[i] == 0 and l[i+1] == 0 and l[i+2] == y:
                changed += fill_except(l, y, i, i+1)
            if l[i] == 0 and l[i+1] == y and l[i+2] == 0:
                changed += fill_except(l, y, i, i+2)
            if l[i] == y and l[i+1] == 0 and l[i+2] == 0:
                changed += fill_except(l, y, i+1, i+2)

    return changed

'''
In case four X are alreay filled, and there is another fully filled line
'L' whose four of X posiiton are same as this line, then we cannot X in
the current line where the fifth X in 'L' exist.

Otherwise, the current line and the line L becomes same, and violates
the rule 4.
'''
def compare(l, L):
    if l.count(0) == 0 or L.count(0) > 0:
        return 0

    changed = 0
    l1, l2 = bucket(l)
    L1, L2 = bucket(L)

    if len(l1) == 4 and l1.issubset(L1):
        for i in L1.difference(l1):
            if l[i] == 0:
                changed += 1
                l[i] = 2

    if len(l2) == 4 and l2.issubset(L2):
        for i in L2.difference(l2):
            if l[i] == 0:
                changed += 1
                l[i] = 1

    return changed

def bucket(l):
    b1 = set()
    b2 = set()
    for i in range(len(l)):
        if l[i] == 1:
            b1.add(i)
        elif l[i] == 2:
            b2.add(i)

    return b1, b2

main()
