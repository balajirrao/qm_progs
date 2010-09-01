#!/usr/bin/python

# Quantum walk simulation

import qwalk_lib as qwalk
import density_matrix as dm
import discord
import sys
import numpy as np
import matplotlib.pyplot as plt

# Graphical output needed ?
graph = True

# Print the probs ?
print_probs = False

step_spec = sys.argv[1].split(':')

step_start = step_end = step_steps = 0

if (len(step_spec) == 0 or len(step_spec) > 3) :
	print "Bad args!"
	exit()
else :
	if len(step_spec) >= 1 :
		step_start = int(step_spec[0])
		step_end = step_start
		step_steps = 1
		
	if len(step_spec) >= 2 :
		step_end = int(step_spec[1])
		step_steps = 10

	if len(step_spec) == 3 :
		step_steps = int(step_spec[2])

	#Noise parameter
	lamda = 0

#Coinspace is invariant
coin_space = qwalk.CoinSpace()
dim_A = len(coin_space)

values = list()
steps_space = range(step_start, step_end, step_steps)
for steps in steps_space:
	pos_list = range (-steps, steps + 1)
	pos_space = qwalk.PositionSpace(pos_list)
	dim_B = len(pos_space)

	E = [np.mat(np.empty((2 * (2 * steps + 1), 2 * (2 * steps + 1)))) for i in range (2)]
	E[0] = np.kron(np.mat([[ np.sqrt(1 - lamda), 0], [0, 1]]), pos_space.I)
	E[1] = np.kron(np.mat([[0, 0], [np.sqrt(lamda), 0]]), pos_space.I)

	walk = qwalk.QWalk(pos_space, coin_space)

	rho_CP = walk.do(steps, E)
	_v = discord.T_min(rho_CP, dim_A, dim_B, qubit_subsys = 1) - dm.S_cond(rho_CP, (2, 1), dim_A, dim_B)
	values.append(_v)

	print steps, " : ", _v

#rho_P = dm.partial_trace(rho_CP, subsys = 1, dim_B = len(pos_space), dim_A = len(coin_space))
#print np.diag(rho_P), np.trace(rho_P)

plt.plot(steps_space, values)
plt.xlabel('Steps')
plt.ylabel('D')
plt.show()
