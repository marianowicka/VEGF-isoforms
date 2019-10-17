import numpy as np

res = np.genfromtxt('bun_end.txt')

n = 181
m = 999


def com_val(ind):
	ret = []
	for j in xrange(0,n):
		parM = []
		i = j
		while i <m*n:
			parM.append(res[i])
			i = i+n
		parM = np.array(parM)
		ret.append(np.percentile(parM, ind, axis=0))
	return(np.array(ret))

x = com_val(5)
per5 = np.array(x)
x = com_val(50)
per50 = np.array(x)
x = com_val(95)
per95 = np.array(x)

np.savetxt('per5_end.txt', per5)
np.savetxt('per50_end.txt', per50)
np.savetxt('per95_end.txt', per95)


#	print 'mean', np.mean(parM, axis = 0)
#	print 'median', np.median(parM, axis = 0)
#	print '5th', np.percentile(parM, 5, axis=0)
#	print '95th', np.percentile(parM, 95, axis=0)
#	print '------------------------------'

