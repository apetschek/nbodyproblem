import numpy as np; import scipy as sc ;from time import time; import simplejson
from solvers import rk4, rk5, rk8, leapfrog2, leapfrog4, rk_adaptive
from helper_functions import accel, masses, positions, velocities, progress
import matplotlib.pyplot as plt

# System parameters
M_sun = 1.98 * 10**30       
AU    = 1.496 * 10**10      
M_param = M_sun        
R_param = AU*.1            
V_param = 100000             
N = 100                
h = 250              
steps = 100*105000/2500 #one orbit is 105,000 steps @ h = 0.1 = 10 sec to run leapfrog.

# Initialize initial conditions
m = masses(N,M_param); r = positions(N,R_param); v = velocities(N,V_param)   
# r[0][0] = .5*R_param; r[1][0] = 0; r[2][0] = 0; r[0][1] = -.5*R_param; r[1][1] = 0; r[2][1] = 0
# v[0][0] = 0; v[1][0] = 125000; v[2][0] = 0; v[0][1] = 0; v[1][1] = -125000; v[2][1] = 0
Rx = [[] for i in range(N)]; Ry = [[] for i in range(N)]; Rz = [[] for i in range(N)] 
Vx = [[] for i in range(N)]; Vy = [[] for i in range(N)]; Vz = [[] for i in range(N)] 

t0 = time()  
for step in range(steps+1):
	progress(step,steps,t0)

	######################################################################
	#                   SELECT DESIRED SOLVER BELOW                      #
	#                                                                    #
	# Runge-Kutta Methods                                                #
	#v,r = rk4(accel,m,r,h,v)             #4th Order (fixed time step)   # 
	#v,r = rk5(accel,m,r,h,v)			  #5th Order (fixed time step)   #
	#v,r = rk8(accel,m,r,h,v)
	#v,r,h = rk_adaptive(accel,m,r,h,v,0) #5th Order (Adaptive time step)#
	#																	 #
	# Symplectic Methods											     #
	#v,r = leapfrog2(accel,m,r,h,v)   #leapfrog 2nd Order (work)s        #
	v,r = leapfrog4(accel,m,r,h,v)   #leapfrog 4th Order  (works)       #
	#																	 #
	# Multistep Methods													 #
	######################################################################	
	
	# Store position
	for i in range(N):
		Rx[i].append(r[0][i])
		Ry[i].append(r[1][i]) 
		Rz[i].append(r[2][i])

	# Store velocity
	for i in range(N):
		Vx[i].append(v[0][i])
		Vy[i].append(v[1][i]) 
		Vz[i].append(v[2][i])

t1 = time()
print t1-t0 

# Save output
with open ('N_output_LF4_250s_BIG.json', 'w') as n:
	simplejson.dump(N, n)
with open ('Rx_output_LF4_250s_BIG.json', 'w') as rx:
	simplejson.dump(Rx, rx)
with open ('Ry_output_LF4_250s_BIG.json', 'w') as ry:
	simplejson.dump(Ry, ry)
with open ('Rz_output_LF4_250s_BIG.json', 'w') as rz:
	simplejson.dump(Rz, rz)
with open ('Vx_output_LF4_250s_BIG.json', 'w') as vx:
	simplejson.dump(Vx, vx)
with open ('Vy_output_LF4_250s_BIG.json', 'w') as vy:
	simplejson.dump(Vy, vy)
with open ('Vz_output_LF4_250s_BIG.json', 'w') as vz:
	simplejson.dump(Vz, vz)
