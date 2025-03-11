# making a neural network trained with the genetic algorithm

import ctypes


class ga_conf(ctypes.Structure):
	_fields_ = [("ranges"), ctypes.pointer(ctypes.c_float),
				("dims"), ctypes.c_int,
				("size"), ctypes.c_int,
				("tour_size"), ctypes.c_int,
				("mut_rate"), ctypes.c_float,
				("elitis"), ctypes.c_float,
				("gens"), ctypes.c_int,
				("find_max"), ctypes.c_int,
				("sel_alg"), ctypes.c_int]


# target function type
c_callback = ctypes.CFUNCTYPE(ctypes.pointer(ctypes.c_float))

fun = ctypes.CDLL("libgenetic.so")

fun.myFunction.argtypes = [ctypes.pointer(ga_conf), c_callback]
fun.ga

