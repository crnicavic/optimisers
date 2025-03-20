# making a neural network trained with the genetic algorithm
import numpy as np
import ctypes

def f(x):
	return (x[0] - 2) ** 2 + x[1] ** 2

class ga_conf(ctypes.Structure):
	_fields_ = [
		("ranges", ctypes.POINTER(ctypes.c_float)),
		("dims", ctypes.c_int),
		("size", ctypes.c_int),
		("tour_size", ctypes.c_int),
		("mut_rate", ctypes.c_float),
		("elitis", ctypes.c_float),
		("gens", ctypes.c_int),
		("find_max", ctypes.c_int),
		("sel_alg", ctypes.c_int)
	]

# target function type
ga_callback_type = ctypes.CFUNCTYPE(ctypes.c_float, ctypes.POINTER(ctypes.c_float))
ga_callback = ga_callback_type(f)

libgenetic = ctypes.cdll.LoadLibrary("./libgenetic.so")

libgenetic.ga.argtypes = [ctypes.POINTER(ga_conf), ga_callback_type]
libgenetic.ga.restype = ctypes.POINTER(ctypes.c_float)

global conf
conf = ga_conf()
ranges = np.array([-5, 5, -5, 5], dtype=np.float32)
conf.ranges = ranges.ctypes.data_as(ctypes.POINTER(ctypes.c_float)) 
conf.dims = 2
conf.size = 30
conf.tour_size = 5
conf.mut_rate = 0.4
conf.elitis = 0.1
conf.gens = 30
conf.find_max = 0
conf.sel_alg = 0

res = libgenetic.ga(conf, ga_callback)
print(res[0], res[1])
print(f(res))
