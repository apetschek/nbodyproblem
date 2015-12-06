import numpy as np
import scipy as sc
from time import time
from solvers import rungekutta, leapfrog, leapfrogFR, rkf45
from helper_functions import accel, masses, positions, velocities, progress
import pylab
import simplejson

# System parameters
M_param = 1           #set all masses to m_sun
R_param = 1           #set max initial position to be 1/10*AU
V_param = 1           #set max initial velocity 10**5 m/s

N = 2                 #number of bodies in system
h = .001              #time step: (1 sec)
steps = 2000         #number of time steps
 
# Initialize initial conditions
m = masses(N,M_param); r = positions(N,R_param); v = velocities(N,V_param)   

r[0][0] = 0.5; r[1][0] = 0; r[2][0] = 0
r[0][1] = 0; r[1][1] = 0; r[2][1] = 0

v[0][0] = 0; v[1][0] = 1; v[2][0] = 0
v[0][1] = 0; v[1][1] = -1; v[2][1] = 0

Rx = [[] for i in range(N)]; Ry = [[] for i in range(N)]; Rz = [[] for i in range(N)] 
Vx = [[] for i in range(N)]; Vy = [[] for i in range(N)]; Vz = [[] for i in range(N)] 
norms = []

t0 = time()  

for step in range(steps):
	progress(step,steps,t0)

	######################################################################
	#                   SELECT DESIRED SOLVER BELOW                      #
	#                                                                    #
	#v,r = rungekutta(accel,m,r,h,v)   #RK 4th Order                     # 
	v,r = leapfrog(accel,m,r,h,v)     #leapfrog 2nd Order                #
	#v,r = leapfrogFR(accel,m,r,h,v)  #leapfrog 4th Order                #
	#v,r,h = rkf45(accel,m,r,h,v,0)   #RK 5th Order w/ Adaptive timestep #
	######################################################################	
	
	# Store position
	for i in range(N):
		Rx[i].append(r[0][i])
		Ry[i].append(r[1][i]) 
		Rz[i].append(r[2][i])

	# Save norms for density as fct of radius calculations
	if step == 0:
		norms_tmp=[]
		for i in range(N):
			norms_tmp.append(np.sqrt((r[0][i])**2 + (r[1][i])**2 + (r[2][i])**2))
		norms.append(norms_tmp)
	
	if (step % 5000 == 0):
			norms_tmp = []
			for n in range(N):
				norms_tmp.append(np.sqrt((r[0][n])**2 + (r[1][n])**2 + (r[2][n])**2))
			norms.append(norms_tmp)
	
	# Store velocity
	for i in range(N):
		Vx[i].append(r[0][i])
		Vy[i].append(r[1][i]) 
		Vz[i].append(r[2][i])

norms_tmp = []
for n in range(N):
	norms_tmp.append(np.sqrt((r[0][n])**2 + (r[1][n])**2 + (r[2][n])**2))
norms.append(norms_tmp)

t1 = time()
print t1-t0                   

# Save output
with open ('testN_output.json', 'w') as n:
	simplejson.dump(N, n)

with open ('testnorms_output.json', 'w') as n:
	simplejson.dump(norms, n)

with open ('testRx_output.json', 'w') as rx:
	simplejson.dump(Rx, rx)

with open ('testRy_output.json', 'w') as ry:
	simplejson.dump(Ry, ry)

with open ('testRz_output.json', 'w') as rz:
	simplejson.dump(Rz, rz)

with open ('testVx_output.json', 'w') as vx:
	simplejson.dump(Vx, vx)

with open ('testVy_output.json', 'w') as vy:
	simplejson.dump(Vy, vy)

with open ('testVz_output.json', 'w') as vz:
	simplejson.dump(Vz, vz)




