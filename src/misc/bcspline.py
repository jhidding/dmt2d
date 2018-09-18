#!/usr/bin/python

from numpy import *
from numpy import linalg
from numpy import random
import pylab


class BC:
	def __init__(self, I):
		self.I = I
		self.Ix = roll(self.I, 1, 0) - self.I
		self.Iy = roll(self.I, 1, 1) - self.I
		self.Ixy = roll(self.Ix, 1, 1) - self.Ix

		self.A = [
			[1, 0, 0, 0,	0, 0, 0, 0,	0, 0, 0, 0,	0, 0, 0, 0],
			[1, 1, 1, 1,	0, 0, 0, 0,	0, 0, 0, 0,	0, 0, 0, 0],
			[1, 0, 0, 0,	1, 0, 0, 0,	1, 0, 0, 0,	1, 0, 0, 0],
			[1, 1, 1, 1,	1, 1, 1, 1,	1, 1, 1, 1,	1, 1, 1, 1],

			[0, 1, 0, 0,	0, 0, 0, 0,	0, 0, 0, 0,	0, 0, 0, 0],
			[0, 1, 2, 3,	0, 0, 0, 0,	0, 0, 0, 0,	0, 0, 0, 0],
			[0, 1, 0, 0,	0, 1, 0, 0,	0, 1, 0, 0,	0, 1, 0, 0],
			[0, 1, 2, 3,	0, 1, 2, 3,	0, 1, 2, 3,	0, 1, 2, 3],

			[0, 0, 0, 0,	1, 0, 0, 0,	0, 0, 0, 0,	0, 0, 0, 0],
			[0, 0, 0, 0,	1, 1, 1, 1,	0, 0, 0, 0,	0, 0, 0, 0],
			[0, 0, 0, 0,	1, 0, 0, 0,	2, 0, 0, 0,	3, 0, 0, 0],
			[0, 0, 0, 0,	1, 1, 1, 1,	2, 2, 2, 2,	3, 3, 3, 3],

			[0, 0, 0, 0,	0, 1, 0, 0,	0, 0, 0, 0,	0, 0, 0, 0],
			[0, 0, 0, 0,	0, 1, 2, 3,	0, 0, 0, 0,	0, 0, 0, 0],
			[0, 0, 0, 0,	0, 1, 0, 0,	0, 2, 0, 0,	0, 3, 0, 0],
			[0, 0, 0, 0,	0, 1, 2, 3,	0, 2, 4, 6,	0, 3, 6, 9]]

		self.Ai = linalg.inv(self.A)
		print self.Ai

	def __call__(self, x, y):
		X = array([0,0,1,1], int) + int(floor(x))
		p = x - floor(x)
		Y = array([0,1,0,1], int) + int(floor(y))
		q = y - floor(y)

		x = r_[self.I[X,Y], self.Ix[X,Y], self.Iy[X,Y], self.Ixy[X,Y]]
		a = dot(self.Ai, x).reshape(4,4)

		r = 0
		for i in range(4):
			for j in range(4):
				r += a[i,j] * p**i * q**j

		return r
		
if __name__ == "__main__":
	I = random.uniform(0.0, 1.0, (5,5))
	B = BC(I)

	J = zeros((40,40))
	for i in range(40):
		for j in range(40):
			J[i,j] = B(i / 10.0, j / 10.0) 

	pylab.subplot(121)
	pylab.imshow(I, interpolation='nearest')
	pylab.subplot(122)
	pylab.imshow(J, interpolation='nearest')

#	pylab.subplot(223)
#	pylab.imshow(Iy, interpolation='nearest')
#	pylab.subplot(224)
#	pylab.imshow(Ixy, interpolation='nearest')
	
	pylab.show()

