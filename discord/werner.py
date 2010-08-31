#!/usr/bin/python


import numpy as np
import density_matrix
import discord
import matplotlib.pyplot as plt

#Entangled
rho_AB = np.mat('0.5 0 0 0.5; 0 0 0 0; 0 0 0 0; 0.5 0 0 0.5')

ZERO = np.mat('1;0')
ONE = np.mat('0;1')

singlet = (np.kron(ZERO, ZERO) - np.kron(ONE, ONE)) / np.sqrt(2)

bell = [None for i in range(3)]
bell[0] = (np.kron(ONE, ZERO) + np.kron(ZERO, ONE)) / np.sqrt(2)
bell[1] = (np.kron(ONE, ZERO) - np.kron(ZERO, ONE)) / np.sqrt(2)
bell[2] = (np.kron(ZERO, ZERO) + np.kron(ONE, ONE)) / np.sqrt(2)

f_range = np.linspace(0, 1, 5)
values = list()

for f in f_range :
	if (f != 0):
		rho = f * singlet * singlet.H
	else:
		rho = np.zeros((4, 4))

	for i in range(3) :
		rho += ((1 - f) / len(bell)) * bell[i] * bell[i].H

	_x = density_matrix.S_cond(rho, cond=(2,1))
	_y = discord.T_min(rho, qubit_subsys = 1)

	print f, _x, _y

	values.append(_y - _x)

plt.plot(f_range, values)
plt.show()
