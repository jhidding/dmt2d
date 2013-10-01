#!/usr/bin/python

from numpy import *
from numpy import linalg
import sys

A = array([[float(i) for i in line.split()]
        for line in file('bspline.txt', 'r')
        if line[0] != '#'])

M = linalg.inv(A)

for row in M:
        for cel in row:
                print "%2i," % (cel),
        print

