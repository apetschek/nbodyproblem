import random
from time import time
import numpy as np

# Creates mass array
def masses(N,mass):
    return np.ones(N)*mass

# Creates random positions within cube
def positions(N,radius):
    initial_positions = np.ones(3*N)
    for i in range(3*N):
        initial_positions[i] *= random.randrange(-radius,radius) + (random.random()-.5)
    return initial_positions.reshape((3,N))
    
# Creates velocity array
def velocities(N,velocity=0):
    rand = np.ones((3*N))*velocity
    for i in range(len(rand)):
        rand[i] *= np.random.random()
    return rand.reshape((3,N))

# Calculates acceleration in x, y, z directions
def accel(m,r):
    G =-6.67*10**-11
    acceleration = np.zeros([3,len(m)])
    for d in range(3):
        for i in range(len(m)): 
            for j in range(len(m)):
                if j != i:
                    acceleration[d][i] += G*m[j]*(r[d][i]-r[d][j])*(np.linalg.norm(r[:,j]-r[:,i],ord=2))**-3
    return acceleration

# Print simulation progress
def progress(step,steps,t0):
    for i in range(1,10):
        if (step == (steps/10)*i):
            print "{}% Complete".format(int(100*((steps/10)*i)/steps))
            print "{} seconds remain".format(int(((time() - t0) / i)*(10-i)))
