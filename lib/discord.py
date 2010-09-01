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

def T(rho_AB, theta, phi, dim_A = 2, dim_B = 2, qubit_subsys = 2) :
	_basis = [None for i in range(2)]
	_x = np.cos (theta)
	_y = np.sin (theta)
	_basis[0] = np.mat([[_x], [np.exp(1j * phi) * _y]])
	_basis[1] = np.mat([[_y], [-np.exp(1j * phi) * _x]])

	t = 0

	for i in range(len(_basis)):
		if (qubit_subsys == 1) :
			_op = np.kron(projector(_basis[i]), np.eye(dim_B))
		else:
			_op = np.kron(np.eye(dim_A), projector(_basis[i]))

		_rho = _op * rho_AB
#		_rho = _op * rho_AB * _op

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



