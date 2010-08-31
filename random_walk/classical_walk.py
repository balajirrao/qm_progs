#!/usr/bin/python

from numpy import array, zeros, nonzero

from matplotlib.pyplot import plot, show

l = 100

a = zeros((l, 2 * l + 1))

a[0][l] = 1;

for i in range(0, l - 1):
	for j in range (1, 2 * l) :
		a[i + 1][j - 1] += 0.25 * a[i][j];
		a[i + 1][j + 1] += 0.75 * a[i][j];

print sum(a[l - 1])

plot(range(0, l), a[l - 1][(nonzero(a[l - 1]))])

show()
