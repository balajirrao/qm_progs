#!/usr/bin/python

# Quantum walk simulation

from __future__ import print_function

import qwalk_lib as qwalk
import density_matrix as dm
import sys
import numpy as np
import matplotlib.pyplot as plt
from optparse import OptionParser

import discord
import qmid

def print_v(*args) :
	global options

	if (options.verbose > 0) :
		print (*args)
	

usage = 'usage %prog : [options]'
parser = OptionParser(usage = usage)

parser.add_option('', '--qd', help = 'Calculate QD', action='store_true', dest = 'qd_needed')
parser.add_option('', '--qmid', help = 'Calculate QMID', action='store_true', dest = 'qmid_needed')

parser.add_option('-N', '--node-range', type = 'int',  nargs = 3, help = 'Range of nodes to use', dest = 'node_range', metavar = "START END RESOLUTION")
parser.add_option('-n', '--nodes', type = 'int',  default = 25, help = 'Number of nodes', dest = 'nodes', metavar = "NODES")

parser.add_option('-L', '--lamda-range', type = 'float',  nargs = 3, help = 'Range of lamdas to use for amp damping', dest = 'lamda_range', metavar = "START END RESOLUTION")
parser.add_option('-l', '--lamda', type = 'float',  default = 0, help = 'Lamda param for amp damping', dest = 'lamda', metavar = "LAMDA")

parser.add_option('-t', '--traversals', type = 'int',  default = 5, help = 'Number of times to traverse the circle', dest = 'traversals', metavar = "N")
parser.add_option('-T', '--traversals-range', type = 'int',  nargs = 3, help = 'Number of times to traverse the circle', dest = 'traversals_range', metavar = "N")

parser.add_option('-g', '--graphical', action = 'store_true', help = 'show graphical output', default = True, dest = 'graphical')
parser.add_option('-v', '--verbose', help = 'verbose output ', type = int, default = 1, dest = 'verbose', metavar = "LEVEL")

(options, args) = parser.parse_args()

#Noise parameter
lamda = 0

#Coinspace is invariant
coin_space = qwalk.CoinSpace()
dim_A = len(coin_space)

output_list_dict = dict()

if (options.qd_needed) :
	output_list_dict['QD'] = list()

if (options.qmid_needed) :
	output_list_dict['QMID'] = list()

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
			print_v ("Nodes : ", nodes, " Lambda = ", lamda, 'Traversals = ', traversals)

			rho_CP = walk.do(traversals * nodes, E)
			if (options.qd_needed) :
				_qd = discord.T_min(rho_CP, dim_A, dim_B, qubit_subsys = 1) - dm.S_cond(rho_CP, (2, 1), dim_A, dim_B).real
				output_list_dict['QD'].append(_qd)
				print_v ("\tQD = ", _qd)

			if (options.qmid_needed) :
				_qm = qmid.qmid(rho_CP, dim_A, dim_B)
				output_list_dict['QMID'].append(_qm)
				print_v ("\tQMID = ", _qm)
			
			if (options.verbose > 1) :
				rho_P = dm.partial_trace(rho_CP, subsys = 1, dim_B = len(pos_space), dim_A = len(coin_space))
				print_v ("\tProbabilities : ", np.diag(rho_P), ". Sum = ", np.trace(rho_P))

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

	for k in output_list_dict.keys() :
		_y = output_list_dict[k]
		plt.plot(_x, _y, label = k)

	plt.legend()
	plt.xlabel(_label)
	plt.title(_title)

	plt.ylabel('Y')
	plt.show()

