#!/usr/bin/python

# Quantum walk simulation

from numpy import mat, zeros, empty, asscalar, nonzero, kron, identity, sqrt, arange, array, pi
from numpy import cos, sin
from matplotlib.pyplot import plot, show, xticks

import sys

def outer(v1, v2) :
	return v1 * v2.H;

def projector(vec) :
	return outer(vec, vec);

# Returns a position vector of appropriate length
def pos_vec(pos) :
	global position_vectors
	
	return position_vectors[pos]

def prob(v, M) :
	return asscalar(v.H * (M * v));

# Graphical output needed ?
graph = True

# Print the probs ?
print_probs = False

if (len(sys.argv) < 3):
	print "Bad Bad!"
	print "Usage : ", sys.argv[0], " circle_n circle_n"
	sys.exit()

circle_n = int(sys.argv[1])
steps = int(sys.argv[2])

# Coin head and tail
C_H = mat('1 ; 0', dtype=int)
C_T = mat('0 ; 1', dtype=int)

# Hadamard
#H = mat('1 1;1 -1', dtype=int) / sqrt(2);

theta = pi / 4

H = mat([[cos (theta), sin (theta)], [sin (theta), - cos(theta)]])

# Identity for coin and position spaces
I_C = identity(2)
I_P = identity(circle_n)

#pre create position vectors
position_vectors = [mat(zeros((circle_n, 1), dtype=int)) for i in range (circle_n)]

# Initialize position vectors
for i in range(len(position_vectors)):
	position_vectors[i][i] = 1

# Operators for moving Left and Right
L_op = zeros((circle_n, circle_n), dtype=int)
for i in range(0, circle_n) :
	L_op += outer(pos_vec(i - 1), pos_vec(i))

R_op = zeros((circle_n, circle_n), dtype=int)
for i in range(-1, circle_n - 1) :
	R_op += outer(pos_vec(i + 1), pos_vec(i))

S_op = kron(projector(C_H), R_op) + kron(projector(C_T), L_op)
Op = S_op * kron(H, I_P)

# List of projectors for position
M = [mat(empty((2 * circle_n, circle_n), dtype=int)) for i in range (circle_n)]
for j in range(0, circle_n) :
	M[j] = (kron(I_C, projector(pos_vec(j))))

# initial state 
#init = kron(C_T, pos_vec(0)) # Asymmetric walk
#init = kron(C_H, pos_vec(0)) # Asymmetric walk
init = kron (C_H + 1j * C_T, pos_vec(0)) / sqrt(2) # Symmetric Walk


for i in range (0, steps) :
	init =  Op * init
	
	if (print_probs) :
		for j in arange(circle_n -1, -1, -1) :
			print prob(init, M[j]),
		print "\n"

if (graph) :
	values = list();

	for j in arange(circle_n - 1, -1, -1) :
			values.append(prob(init, M[j]))
			

	print 'Sum of probabilities is', abs(sum(values))

	plot(arange(circle_n -1, -1, -1), values)

	# Location of ticks
	loc = range (0, circle_n, circle_n / 5)
	xticks(loc)
	
	show()
