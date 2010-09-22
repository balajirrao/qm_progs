#!/usr/bin/python

#Density Matrix Library

import numpy as np
import itertools as it

# Index of |ij> <kl|
def rho_idx(i, j, k, l, dim_A = 2, dim_B = 2) :
	return (i * dim_B + j, k * dim_B + l)

def partial_trace(rho_AB, subsys = 2, dim_A = 2, dim_B = 2) :
	d = dim_A * dim_B
	if (rho_AB.shape != (d, d)) :
		print "Bad Matrix"
		raise Exception('Bad Matrix')
	
	if (subsys == 2) :
		dim_reduced = dim_A
	else:
		dim_reduced = dim_B

	rho_reduced = np.zeros((dim_reduced, dim_reduced))

#	idxes = it.product(range(dim_A), repeat=2) # New in 2.6
	idxes = ((x,y) for x in range(dim_reduced) for y in range(dim_reduced))

	for idx in idxes :
		for i in range(d / dim_reduced) :
			x, y = idx
			if (subsys == 1):
				rho_reduced[idx] += rho_AB[rho_idx(i, x, i, y, dim_A, dim_B)]
			else:
				rho_reduced[idx] += rho_AB[rho_idx(x, i, y, i, dim_A, dim_B)]
	
	return np.asmatrix(rho_reduced)

def log2(x):
	if (x == 0):
		return 0
	else :
		return np.log2(x)

def H(ev):
	return -ev * log2(ev)

def S(rho) :
	s = 0
	evals = np.round(np.linalg.eigvals(rho), 5)
	for ev in evals :
		s += H(ev)
	return s

def S_cond(rho_AB, cond = (1, 2), dim_A = 2, dim_B = 2) :
	_x = S(rho_AB) 
	_y = S(partial_trace(rho_AB, cond[0], dim_A, dim_B))

	return _x - _y

def I(rho, dim_A = 2, dim_B = 2) :
	rho_A = partial_trace(rho, subsys = 2, dim_A = dim_A, dim_B = dim_B)
	rho_B = partial_trace(rho, subsys = 2, dim_A = dim_A, dim_B = dim_B)

	return S(rho_A) + S(rho_B) - S(rho)
