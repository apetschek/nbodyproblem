import random
from time import time
import numpy as np

#creates mass array
def masses(N,mass):
    return np.ones(N)*mass

#creates random positions within cube
def positions(N,radius):
    rand = np.ones((3*N))*radius
    for i in range(len(rand)):
        rand[i] *= np.random.random()
    return rand.reshape((3,N))

#creates velocity array
def velocities(N,velocity=0):
    rand = np.ones((3*N))*velocity
    for i in range(len(rand)):
        rand[i] *= np.random.random()
    return rand.reshape((3,N))

#calculates acceleration in x, y, z directions
def accel(m,r):
    G =-6.67*10**-11
    acceleration = np.zeros([3,len(m)])
    for d in range(3):
        for i in range(len(m)): 
            for j in range(len(m)):
                if j != i:
                    acceleration[d][i] += G*m[j]*(r[d][i]-r[d][j])*(np.linalg.norm(r[:,j]-r[:,i],ord=2))**-3
    return acceleration

def progress(step,steps,t0):
    tenpercent = steps/10
    for i in range(1,10):
        if (step == tenpercent*i):
            print "{}% Complete".format(int(tenpercent*i*.01))
            print "{} seconds remain".format(int(((time() - t0) / i)*(10-i)))
