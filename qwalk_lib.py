#Module for Quantum Walks
import numpy as np

#Utility functions
def projector(vec) :
	return np.outer(vec, vec);

tensor_product = np.kron

class PositionSpace(list) :
	__slots__ = 'Rop', 'Lop', 'I', '__vectors', '__pos_list'

	def __new__(cls, *args, **kwrds):
		return list.__new__(cls)

	def __init__(self, pos_list) :
		list.__init__(self)

		# Create standard basis vectors for position space
		m = np.asmatrix(np.eye(len(pos_list), dtype=int))
		for bra in m:
			list.append(self, bra.H)

		self.__pos_list = pos_list

		#Identity operator
		self.I = np.eye(len(self.__pos_list))

	def __getitem__(self, index) :
		_idx = self.__pos_list.index(index)
		return list.__getitem__(self, _idx);

	def __index__(self, item) :
		print item
		return list.index(self, item)

class CoinSpace(dict) :
	__slots__ = 'I'

	def __new__(cls):
		return dict.__new__(cls)

	def __init__(self):
		_head = np.mat('1 ; 0')
		_tail = np.mat('0 ; 1')

		dict.__init__(self, {'H' : _head, 'T' : _tail})

		self.I = np.eye(2)

	def flip_op(self, theta = np.pi / 4):
		return np.mat([[np.cos (theta), np.sin (theta)], [np.sin (theta), - np.cos(theta)]])


class QWalk:
	def __init__(self, position_space, coin_space, circle=False, symmetric_init=True) :

		self.__pos_space = position_space
		self.__coin_space = coin_space
		self.__circle = circle
		self.__symmetric_init = symmetric_init

		# Create the Right and Left operators for this space
		self.__Rop = np.asmatrix(np.eye(len(position_space), k=-1))
		self.__Lop = np.asmatrix(np.eye(len(position_space), k=+1))

		if (circle == True) :
			self.__Rop[0][-1] = 1
			self.__Lop[-1][0] = 1
		
		_move_op = \
			tensor_product(projector(self.__coin_space['H']), self.__Rop) + \
			tensor_product(projector(self.__coin_space['T']), self.__Lop)

		self.__walk_op = \
			_move_op * tensor_product(self.__coin_space.flip_op(), self.__pos_space.I)
		
		if (self.__symmetric_init == False) :
			self.__rho = projector(tensor_product(self.__coin_space['T'], self.__pos_space[0]))
		else :
			self.__rho = projector(tensor_product(
				self.__coin_space['T'] + 1j * self.__coin_space['H'],
				self.__pos_space[0]))

	def do(self, steps, operation_map=[]):
		if (steps <= 0) :
			print "Steps should be greater than zero"
			return


		for _i in range(steps) :
			self.__rho = (self.__walk_op * self.__rho) * self.__walk_op.H;
			_rho = np.zeros_like(self.__rho)

			for E in operation_map:
				_rho += (E * self.__rho) * E.H

			if (len(operation_map) != 0) :
				self.__rho = _rho
		
		return self.__rho
