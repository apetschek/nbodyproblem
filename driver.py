import numpy as np; import scipy as sc ;from time import time; import simplejson
from solvers import rk4, rk5, rk8, leapfrog2, leapfrog4, rk_adaptive
from helper_functions import accel, masses, positions, velocities, progress
import matplotlib.pyplot as plt

# System parameters
M_sun = 1.98 * 10**30       
AU    = 1.496 * 10**10      
M_param = M_sun        
R_param = AU            
V_param = 100000             
N = 6                
h = 250             
steps = 10000

# Initialize initial conditions
m = masses(N,M_param); r = positions(N,R_param); v = velocities(N,V_param)   

Rx = [[] for i in range(N)]; Ry = [[] for i in range(N)]; Rz = [[] for i in range(N)] 
Vx = [[] for i in range(N)]; Vy = [[] for i in range(N)]; Vz = [[] for i in range(N)] 

t0 = time()  
for step in range(steps):
	progress(step,steps,t0, N, Rx, Ry, Rz, Vx, Vy, Vz)

	######################################################################
	#                   SELECT DESIRED SOLVER BELOW                      #
	#                                                                    #
	# Runge-Kutta Methods                                                #
	#v,r = rk4(accel,m,r,h,v)             #4th Order (fixed time step)   # 
	#v,r = rk5(accel,m,r,h,v)			  #5th Order (fixed time step)   #
	#v,r = rk8(accel,m,r,h,v) 											 #
	#v,r,h = rk_adaptive(accel,m,r,h,v,0) #5th Order (Adaptive time step)#
	#																	 #
	# Symplectic Methods											     #
	v,r = leapfrog2(accel,m,r,h,v)    #2nd Order (fixed time step        #
	#v,r = leapfrog4(accel,m,r,h,v)   #4th Order (fixed time step)       #
	#																	 #
	# 																	 #
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
with open ('N_output.json', 'w') as n:
	simplejson.dump(N, n)
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
