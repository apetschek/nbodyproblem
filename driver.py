import numpy as np
import scipy as sc
from time import time
from solvers import rungekutta, leapfrog, leapfrogFR, rkf45
from helper_functions import accel, masses, positions, velocities, progress
import pylab
import simplejson

# System parameters
M_sun = 1.98 * 10**30       #mass of sun in kg
AU    = 1.496 * 10**10      #distance from earth to sun in m
M_param = M_sun             #set all masses to m_sun
R_param = AU * .1           #set max initial position to be 1/10*AU
V_param = 100000            #set max initial velocity 10**5 m/s

N = 30                      #number of bodies in system
h = 1                       #time step: (1 sec)
steps = 20000               #number of time steps
 
# Initialize initial conditions
m = masses(N,M_param); r = positions(N,R_param); v = velocities(N,V_param)   
Rx = [[] for i in range(N)]; Ry = [[] for i in range(N)]; Rz = [[] for i in range(N)] 
Vx = [[] for i in range(N)]; Vy = [[] for i in range(N)]; Vz = [[] for i in range(N)] 
norms = [[] for i in range(5)]

t0 = time()  

for step in range(steps):
	progress(step,steps,t0)

	######################################################################
	#                   SELECT DESIRED SOLVER BELOW                      #
	#                                                                    #
	#v,r = rungekutta(accel,m,r,h,v)  #RK 4th Order                      # 
	v,r = leapfrog(accel,m,r,h,v)     #leapfrog 2nd Order                #
	#v,r = leapfrogFR(accel,m,r,h,v)  #leapfrog 4th Order                #
	#v,r,h = rkf45(accel,m,r,h,v,0)   #RK 5th Order w/ Adaptive timestep #
	######################################################################	
	
	# Store position
	for i in range(N):
		Rx[i].append(r[0][i])
		Ry[i].append(r[1][i]) 
		Rz[i].append(r[2][i])

	if step == 0:
		for i in range(N):
			norms[0].append(np.sqrt((r[0][i])**2 + (r[1][i])**2 + (r[2][i])**2))
	
	for i in range(4):
		if (step == (steps/4)*i):
			for n in range(N):
				norms[i].append(np.sqrt((r[0][n])**2 + (r[1][n])**2 + (r[2][n])**2))

	# Store velocity
	for i in range(N):
		Vx[i].append(r[0][i])
		Vy[i].append(r[1][i]) 
		Vz[i].append(r[2][i])


for n in range(N):
	norms[4].append(np.sqrt((r[0][n])**2 + (r[1][n])**2 + (r[2][n])**2))

t1 = time()
print t1-t0                   

# Save output
with open ('N_output.json', 'w') as n:
	simplejson.dump(N, n)

with open ('norms_output.json', 'w') as n:
	simplejson.dump(norms, d)

with open ('Rx_output.json', 'w') as rx:
	simplejson.dump(Rx, rx)

with open ('Ry_output.json', 'w') as ry:
	simplejson.dump(Ry, ry)

with open ('Rz_output.json', 'w') as rz:
	simplejson.dump(Rz, rz)

with open ('Vx_output.json', 'w') as vx:
	simplejson.dump(Vx, vx)

with open ('Vy_output.json', 'w') as vy:
	simplejson.dump(Vy, vy)

with open ('Vz_output.json', 'w') as vz:
	simplejson.dump(Vz, vz)




