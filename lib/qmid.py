import numpy as np
import density_matrix
import matplotlib.pyplot as plt

def proj(v) :
	return v * v.H

def qmid(rho, dim_A = 2, dim_B = 2) :
	rho_A = density_matrix.partial_trace(rho, 2, dim_A, dim_B)
	rho_B = density_matrix.partial_trace(rho, 1, dim_A, dim_B)

	val_A, vec_A = np.linalg.eig(rho_A)
	val_B, vec_B = np.linalg.eig(rho_B)

	pi = list()
	for v_A in vec_A :
		for v_B in vec_B:
			pi.append(np.kron(proj(v_A.T), proj(v_B.T)))

	_rho = np.zeros_like(rho)
	for p in pi :
		_rho += p * rho * p

	return density_matrix.I(rho, dim_A, dim_B) - density_matrix.I(_rho, dim_A, dim_B)
