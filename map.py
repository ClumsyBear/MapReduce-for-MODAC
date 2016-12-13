#!/sw/lsa/centos7/python-anaconda2/201607/bin/python

# #!/usr/bin/env python2.7 

import fileinput
import numpy as np
from numpy.linalg import inv
from sklearn import linear_model

X = []
y = []
sid = []
count = 0

for line in fileinput.input():
	rowitems = line.rstrip('\n').split('\t')
	X_item = map(float, rowitems[0:-2])
	y_item = float(rowitems[-2])
	sid_item = int(rowitems[-1])
	X.append(X_item)
	y.append(y_item)
	sid.append(sid_item)
	count += 1

# print count
X_train = np.array(X)
y_train = np.array(y)
sid_train = np.array(sid)

# regr = linear_model.LinearRegression()            # this runs linear regression
regr = linear_model.LassoLarsIC(criterion='bic', fit_intercept=False)  # this runs lasso regression

regr.fit(X_train, y_train)

fit_coef = regr.coef_
# P   # for GLM, P is diagonal variance matrix estimated using fit_coef_
fit_covinv = np.dot(X_train.T, X_train)  # for GLM, use np.dot(np.dot(X_train.T, P), X_train) instead
fit_lambda = regr.alpha_

A = np.dot(inv(fit_covinv), X_train.T)
e = y - np.dot(X_train, fit_coef)  # for GLM, transform np.dot(X_train, fit_coef) using link function
betahat_c = fit_coef + np.dot(A, e)

local_id = np.random.randint(1e6)  # sub-dataset stored in different machines have different id

print str(local_id) + '\t' + 'b' + '\t' + '\t'.join(map(str, betahat_c))
row_count = 0
for row in fit_covinv:
	print str(local_id) + '\t' + str(row_count) + '\t' + '\t'.join(map(str, row))
	row_count += 1
print str(local_id) + '\t' + 'c' + '\t' + str(count)
print str(local_id) + '\t' + 'l' + '\t' + str(fit_lambda)



