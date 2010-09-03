#!/usr/bin/python

# Quantum walk simulation

import qwalk_lib as qwalk
import density_matrix as dm
import discord
import sys
import numpy as np
import matplotlib.pyplot as plt
from optparse import OptionParser

usage = 'usage %prog : [options]'
parser = OptionParser(usage = usage)

parser.add_option('-g', '--graphical', action = 'store_true', help = 'show graphical output', default = True, dest = 'graphical')
parser.add_option('-v', '--verbose', help = 'verbose output ', type = int, default = 1, dest = 'verbose', metavar = "LEVEL")
parser.add_option('-N', '--node-range', type = 'int',  nargs = 3, help = 'Range of nodes to use', dest = 'node_range', metavar = "START END RESOLUTION")
parser.add_option('-n', '--nodes', type = 'int',  help = 'Number of nodes', dest = 'nodes', metavar = "NODES")
parser.add_option('-t', '--traversals', type = 'int',  default = 1, help = 'Number of times to traverse the circle', dest = 'traversals', metavar = "N")

(options, args) = parser.parse_args()

#Noise parameter
lamda = 0

#Coinspace is invariant
coin_space = qwalk.CoinSpace()
dim_A = len(coin_space)

values = list()

nodes_space = range(*options.node_range)

for nodes in nodes_space:
	pos_list = range (-nodes, nodes + 1)
	pos_space = qwalk.PositionSpace(pos_list)
	dim_B = len(pos_space)

	E = [np.mat(np.empty((2 * (2 * nodes + 1), 2 * (2 * nodes + 1)))) for i in range (2)]
	E[0] = np.kron(np.mat([[ np.sqrt(1 - lamda), 0], [0, 1]]), pos_space.I)
	E[1] = np.kron(np.mat([[0, 0], [np.sqrt(lamda), 0]]), pos_space.I)

	walk = qwalk.QWalk(pos_space, coin_space, circle = True)

	rho_CP = walk.do(options.traversals * nodes)

	_qd = discord.T_min(rho_CP, dim_A, dim_B, qubit_subsys = 1) - dm.S_cond(rho_CP, (2, 1), dim_A, dim_B).real
	values.append(_qd)

	if options.verbose > 0 :
		print "Nodes : ", nodes
		print "\tQD = ", _qd
		
		if options.verbose > 1 :
			rho_P = dm.partial_trace(rho_CP, subsys = 1, dim_B = len(pos_space), dim_A = len(coin_space))
			print "\tProbabilities : ", np.diag(rho_P), ". Sum = ", np.trace(rho_P)

if (options.graphical) :
	plt.plot(nodes_space, values)
	plt.xlabel('Steps')
	plt.ylabel('D')
	plt.show()
