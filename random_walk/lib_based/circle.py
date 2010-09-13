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
parser.add_option('-n', '--nodes', type = 'int',  default = 25, help = 'Number of nodes', dest = 'nodes', metavar = "NODES")

parser.add_option('-L', '--lamda-range', type = 'float',  nargs = 3, help = 'Range of lamdas to use for amp damping', dest = 'lamda_range', metavar = "START END RESOLUTION")
parser.add_option('-l', '--lamda', type = 'float',  default = 0, help = 'Lamda param for amp damping', dest = 'lamda', metavar = "LAMDA")

parser.add_option('-t', '--traversals', type = 'int',  default = 5, help = 'Number of times to traverse the circle', dest = 'traversals', metavar = "N")
parser.add_option('-T', '--traversals-range', type = 'int',  nargs = 3, help = 'Number of times to traverse the circle', dest = 'traversals_range', metavar = "N")

(options, args) = parser.parse_args()

#Noise parameter
lamda = 0

#Coinspace is invariant
coin_space = qwalk.CoinSpace()
dim_A = len(coin_space)

values = list()

_options = 0

nodes_space = lamda_space = traversals_space = None

if options.node_range is not None :
	nodes_space = range(*options.node_range)
	_options += 1
else:
	nodes_space = [options.nodes]


if options.lamda_range is not None :
	lamda_space = range(*options.lamda_range)
	_options += 1
else:
	lamda_space = [options.lamda]

if options.traversals_range is not None :
	traversals_space = range(*options.traversals_range)
	_options += 1
else:
	traversals_space = [options.traversals]

if _options == 0 or _options != 1 :
	parser.print_help()
	sys.exit(1)

for nodes in nodes_space:
	if (nodes % 2 == 0) :
		pos_list = range (-nodes / 2, nodes / 2)
	else:
		pos_list = range (-nodes / 2, nodes / 2 + 1)

	pos_space = qwalk.PositionSpace(pos_list)
	dim_B = len(pos_space)

	walk = qwalk.QWalk(pos_space, coin_space, circle = True)

	for lamda in lamda_space:
		E = [np.mat(np.empty((2 * (2 * nodes + 1), 2 * (2 * nodes + 1)))) for i in range (2)]
		E[0] = np.kron(np.mat([[ np.sqrt(1 - lamda), 0], [0, 1]]), pos_space.I)
		E[1] = np.kron(np.mat([[0, 0], [np.sqrt(lamda), 0]]), pos_space.I)

		for traversals in traversals_space :
			rho_CP = walk.do(traversals * nodes, E)
			_qd = discord.T_min(rho_CP, dim_A, dim_B, qubit_subsys = 1) - dm.S_cond(rho_CP, (2, 1), dim_A, dim_B).real
			values.append(_qd)

			if options.verbose > 0 :
				print "Nodes : ", nodes, " Lambda = ", lamda, 'Traversals = ', traversals
				print "\tQD = ", _qd
			
				if options.verbose > 1 :
					rho_P = dm.partial_trace(rho_CP, subsys = 1, dim_B = len(pos_space), dim_A = len(coin_space))
					print "\tProbabilities : ", np.diag(rho_P), ". Sum = ", np.trace(rho_P)

if (options.graphical) :
	if options.node_range is not None :
		_x = nodes_space
		_label = 'Nodes'
		_title = 'Walk on circle with lamda = ' + str(options.lamda) + ', traversals = ' + str(options.traversals)
	elif options.lamda_range is not None :
		_x = lamda_space
		_label = 'Noise Parameter'
		_title = 'Walk on circle with nodes = '+ str(options.nodes) + ', traversals = ' + str(options.traversals)
	else :
		_x = traversals_space
		_label = 'Traversals'
		_title = 'Walk on circle with nodes = ' + str(options.nodes) + ', lamda = ' + str(options.lamda)

	plt.plot(_x, values)
	plt.xlabel(_label)
	plt.title(_title)

	plt.ylabel('D')
	plt.show()

