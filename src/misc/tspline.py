#!/usr/bin/python

from numpy import *
from numpy import linalg
import sys

A = array([[float(i) for i in line.split()]
        for line in file('tspline.txt', 'r')
        if line[0] != '#'])

M = linalg.inv(A)

for row in M:
        for cel in row[:32]:
                print "%3i," % (cel),
        print
        for cel in row[32:]:
                print "%3i," % (cel),
        print

