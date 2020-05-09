from numpy import array, zeros
import matplotlib.pyplot as plt
from timeit import default_timer as timer
from numba import jit


def readarrays(filename):
	#start = timer()
	values = open(filename, "r")
	#print values.read()
	lines = values.readlines()
	#end = timer()
	#print (end-start)


	#Counting
	C = 0
	D = 0
	Dims = []
	A = []

	#start = timer()
	for i in lines:
		if i != "\n":
			D += 1
		if i == "\n":
			C += 1
			Dims.append(D)
			A.append(zeros(D))
			D = 0
	#end = timer()
	#print (end-start)

	#start = timer()
	#Filling
	F = 0
	G = 0
	for i in lines:
		if i != "\n":
			A[F][G] = i
			G += 1
		if i == "\n":
			F += 1
			G = 0
	#end = timer()
	#print (end-start)
	values.close()
	return A,len(A)

def readmatrices(filename):
	#start = timer()
	values = open(filename, "r")
	#print values.read()
	lines = values.readlines()
	#end = timer()
	#print (end-start)


	#Counting
	C = 0
	D = 0
	Dims = []
	A = []

	#start = timer()

	for i in range(len(lines)):
		if lines[i] != "\n":
			D += 1
		if lines[i] == "\n":
			C += 1
			Dims.append(D)
			A.append(zeros(shape=(D,len(lines[i-1].split()))))
			D = 0


	#start = timer()
	#Filling
	F = 0
	G = 0
	for i in lines:
		if i != "\n":
			for j in range(len(i.split())):
				A[F][G][j] = i.split()[j]
			G += 1
		if i == "\n":
			F += 1
			G = 0
	#end = timer()
	#print (end-start)
	values.close()

	return A,len(A)
