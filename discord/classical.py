#!/usr/bin/python


import numpy as np
import density_matrix
import discord
import matplotlib.pyplot as plt

def proj(v) :
	return v * v.H

#Entangled
rho_AB = np.mat('0.5 0 0 0.5; 0 0 0 0; 0 0 0 0; 0.5 0 0 0.5')

ZERO = np.mat('1;0')
ONE = np.mat('0;1')
PLUS = np.mat('1; 1') / np.sqrt(2)
MINUS = np.mat('1; -1') / np.sqrt(2)

psi = np.kron(proj(PLUS), proj (ZERO)) + \
		np.kron(proj(MINUS), proj (PLUS))
#		np.kron(proj(ZERO), proj (MINUS)) + \
#		np.kron(proj(ONE), proj (PLUS))
		
psi /= 2

_x = density_matrix.S_cond(psi)
_y = discord.T_min(psi)

print _y - _x
