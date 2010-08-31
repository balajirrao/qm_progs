#!/usr/bin/python

# Quantum walk simulation

import qwalk_lib as qwalk
import density_matrix as dm
import discord
import sys
import numpy as np

# Graphical output needed ?
graph = True

# Print the probs ?
print_probs = False

steps = int(sys.argv[1])

# Parameter for the Kraus operators
lamda = 0.15

#Inludes 3..

pos_list = range (-steps, steps + 1)
pos_space = qwalk.PositionSpace(pos_list)
coin_space = qwalk.CoinSpace()

E = [np.mat(np.empty((2 * (2 * steps + 1), 2 * (2 * steps + 1)))) for i in range (2)]
E[0] = np.kron(np.mat([[ np.sqrt(1 - lamda), 0], [0, 1]]), pos_space.I)
E[1] = np.kron(np.mat([[0, 0], [np.sqrt(lamda), 0]]), pos_space.I)

walk = qwalk.QWalk(pos_space, coin_space)
rho_CP = walk.do(steps, E)

#rho_P = dm.partial_trace(rho_CP, subsys = 1, dim_B = len(pos_space), dim_A = len(coin_space))
#print np.diag(rho_P), np.trace(rho_P)

dim_A = len(coin_space)
dim_B = len(pos_space)

print discord.T_min(rho_CP, dim_A, dim_B, qubit_subsys = 1) - dm.S_cond(rho_CP, (1, 2), dim_A, dim_B)
