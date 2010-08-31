#!/usr/bin/python


import numpy as np
import density_matrix
import discord
import matplotlib.pyplot as plt

def projector(v) :
	return v * v.H

#Entangled
rho_AB = np.mat('0.5 0 0 0.5; 0 0 0 0; 0 0 0 0; 0.5 0 0 0.5')

ZERO = np.mat('1;0')
ONE = np.mat('0;1')

singlet = projector((np.kron(ZERO, ZERO) - np.kron(ONE, ONE)) / np.sqrt(2))

_x = density_matrix.S_cond(singlet)
_y = discord.T_min(singlet)

print _y - _x
