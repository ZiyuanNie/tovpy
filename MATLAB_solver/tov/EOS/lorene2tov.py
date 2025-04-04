#!/usr/bin/env python

import sys
import fileinput

nlines = -1
for line in fileinput.input():
    if line[0] == '#': continue
    if nlines < 0:
        nlines = int(line.rstrip().split(' ')[0])
        print("#\n"*3,end='')
        print(nlines)
        print("#\n"*3,end='')
        continue
    i, nb, rho, p = [float(x) for x in line.split()]
    i = int(i)
    nb *= 10**(39)
    print("{:d} {:1.7e} {:1.7e} {:1.7e}".format(i, nb, rho, p))
