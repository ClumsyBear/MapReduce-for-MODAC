#!/sw/lsa/centos7/python-anaconda2/201607/bin/python

# #!/usr/bin/env python2.7

import fileinput
import numpy as np
from numpy.linalg import inv

beta_dict = {}
var_dict = {}
size_dict = {}
lambda_dict = {}
p = 0
N = 0.0

for line in fileinput.input():
	# print line.rstrip('\n')  # prints out anything passed to the reducer
	rowitems = line.rstrip('\n').split('\t')
	if rowitems[1] == 'b':
		p = len(rowitems[2:])
		beta_dict[rowitems[0]] = map(float, rowitems[2:])
	elif rowitems[1] == 'c':
		N += int(rowitems[2])
		size_dict[rowitems[0]] = int(rowitems[2])
	elif rowitems[1] == 'l':
		lambda_dict[rowitems[0]] = float(rowitems[2])
	else:
		if rowitems[0] in var_dict.keys():
			pass
			var_dict[rowitems[0]].append(map(float, rowitems[1:]))
		else:
			var_dict[rowitems[0]] = [map(float, rowitems[1:])]

sigbeta_sum = np.array([0.] * p)
sigma_sum = np.zeros((p,p))

for subid in beta_dict.keys():
	beta_dict[subid] = np.array(beta_dict[subid])
	# print beta_dict[subid]
	var_dict[subid].sort(key=lambda x: x[0])
	var_dict[subid] = np.array([item[1:] for item in var_dict[subid]])

for subid in beta_dict.keys():
	sigbeta_sum = sigbeta_sum + np.dot(beta_dict[subid], var_dict[subid])
	sigma_sum = sigma_sum + var_dict[subid]

# print beta_dict
# print var_dict
# print size_dict
# print lambda_dict
# print N

beta_hat = np.dot(inv(sigma_sum), sigbeta_sum)
print beta_hat
print inv(sigma_sum)

