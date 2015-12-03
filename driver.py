import numpy as np
import scipy as sc
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from time import time
from functions import rungekutta, accel, masses, positions, velocities, leapfrog, leapfrogFR, rkf45
import matplotlib.animation as animation
import pylab

######################
# INITIAL CONDITIONS #
######################


N = 3                      #number of bodies in system

M_sun = 1.98*10**30         #mass of sun in kg
AU = 1.496*10**10           #distance from earth to sun in m
M_param = M_sun             #set all masses to m_sun
R_param = AU*.01            #set max initial position to be 1/10*AU
V_param = 100000            #set max initial velocity 10**5 m/s

m = masses(N,M_param)       #create masses vector for N bodies
r = positions(N,R_param)    #create positions vector 
v = velocities(N,V_param)   #create velocities vector 

Rx = [[] for i in range(N)] #create list of x positions to be filled in for plotting
Ry = [[] for i in range(N)] #create list of y positions
Rz = [[] for i in range(N)] #create list of z positions

h = 1                       #time step: (1 sec)
steps = 100                #number of time steps

#########
# SOLVE #
#########

t0 = time()                 #initialize timer

for i in range(steps):

	######################################################################
	#                   SELECT DESIRED SOLVER BELOW                      #
	#                                                                    #
	#v,r = rungekutta(accel,m,r,h,v)  #RK 4th Order                      # 
	v,r = leapfrog(accel,m,r,h,v)     #leapfrog 2nd Order                #
	#v,r = leapfrogFR(accel,m,r,h,v)  #leapfrog 4th Order                #
	#v,r,h = rkf45(accel,m,r,h,v,0)   #RK 5th Order w/ Adaptive timestep #
	######################################################################	
	
	for i in range(N):
		Rx[i].append(r[0][i]) #store x values per body
		Ry[i].append(r[1][i]) #store y values
		Rz[i].append(r[2][i]) #store z values

t1 = time()
print t1-t0                   #print computation time

########
# PLOT #
######## 


# fig = plt.figure(figsize=(10,10))
# ax = fig.add_subplot(111, projection = '3d')
# for i in range(N):
# 	ax.plot(Rx[i], Ry[i], Rz[i])
# ax.set_xlim(-R_param,R_param)
# ax.set_ylim(0,R_param)
# ax.set_zlim(0,R_param)
# ax.set_xlabel('X')
# ax.set_ylabel('Y')
# ax.set_zlabel('Z')
# plt.show()


fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(111, projection = '3d')
for j in range(steps):
	for i in range(N):
		pylab.scatter(Rx[i][j], Ry[i][j], Rz[i][j])
		plt.draw()
ax.set_xlim(-R_param,R_param)
ax.set_ylim(0,R_param)
ax.set_zlim(0,R_param)
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')





