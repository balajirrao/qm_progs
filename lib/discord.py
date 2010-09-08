#!/usr/bin/python

import numpy as np
import density_matrix
import matplotlib.pyplot as plt
import multiprocessing as mp

res = 25
thread_count = 2

def projector(v) :
	return v * v.H

def _T(args) :
	return T(*args)

_X = np.mat('0 1; 1 0')
_Y = np.mat('0 -1j; 1j 0')
_Z = np.mat('1 0; 0 -1')
_I = np.eye(2)

def E(n) :
	_nx, _ny, _nz = n

	return 0.5 * (_I + _nx * _X + _ny * _Y + _nz * _Z)

def T(rho_AB, theta, phi, dim_A = 2, dim_B = 2, qubit_subsys = 2) :

	_n = [None for i in range(3)]
	_n[0] = (np.sin (theta) * np.cos (phi), \
			np.sin (theta) * np.sin (phi), \
			np.cos (theta))
	
	_n[1] = (np.sin (theta - np.pi) * np.cos (phi), \
			np.sin (theta - np.pi) * np.sin (phi), \
			np.cos (theta - np.pi))

	_n[2] = (np.sin (np.pi - theta) * np.cos (np.pi + phi), \
			np.sin (np.pi - theta) * np.sin (np.pi + phi), \
			np.cos (np.pi - theta))

	F = [None for i in range(3)]

	F[0] = E(_n[0]);
	F[1] = E(_n[1]) / 2;
	F[2] = E(_n[2]) / 2;

#	print sum(F)

	t = 0

	for i in range(len(F)):
		if (qubit_subsys == 1) :
			_op = np.kron(F[i], np.eye(dim_B))
		else:
			_op = np.kron(np.eye(dim_A), F[i])

		_rho = _op * rho_AB

		p = np.trace(_rho).real
		if (p != 0) :
			_rho_reduced = density_matrix.partial_trace(_rho, qubit_subsys , dim_A = dim_A, dim_B = dim_B) / p
			t += p * density_matrix.S(_rho_reduced)

	return t

def T_min(rho_AB, dim_A = 2, dim_B = 2, qubit_subsys = 2) :
	tpl = ((rho_AB, theta, phi, dim_A, dim_B, qubit_subsys) for theta in np.linspace(0, np.pi, res) for phi in np.linspace(0, 2 * np.pi, res))

#Multiple Threaded
	pool = mp.Pool(thread_count)
	values = pool.map(_T, tpl)

#Single Threaded
#	values = map(_T, tpl)

	return min(values)



