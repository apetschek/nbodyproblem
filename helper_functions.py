import random
from time import time
import numpy as np
import simplejson

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

def energy(m,steps):
    with open('N_output.json', 'r') as n:
        N = simplejson.load(n)
    with open('Rx_output.json', 'r') as rx:
        Rx = simplejson.load(rx)
    with open('Ry_output.json', 'r') as ry:
        Ry = simplejson.load(ry)
    with open('Rz_output.json', 'r') as rz:
        Rz = simplejson.load(rz)
    with open('Vx_output.json', 'r') as vx:
        Vx = simplejson.load(vx)
    with open('Vy_output.json', 'r') as vy:
        Vy = simplejson.load(vy)
    with open('Vz_output.json', 'r') as vz:
        Vz = simplejson.load(vz)

    G=-6.674*10**(-11)
    
    energy_array = []
    
    for step in range(steps):
        system_energy_i = 0
        for i in range(N):
            vx = Vx[i][step]
            vy = Vy[i][step]
            vz = Vz[i][step]
            v = np.sqrt(vx**2 + vy**2 + vz**2)
            ke = 0.5*m[i]*v**2
            pe = 0
            for j in range(N):
                if i!=j:
                    pe += (G*m[i]*m[j])/np.sqrt((Rx[i][step]-Rx[j][step])**2+(Ry[i][step]-Ry[j][step])**2+(Rz[i][step]-Rz[j][step])**2)
            energy=ke+pe
            system_energy_i += energy
        
        energy_array.append(system_energy_i)
    
    return energy_array

def momentum(m,steps): 

    with open('N_output.json', 'r') as n:
        N = simplejson.load(n)
    with open('Rx_output.json', 'r') as rx:
        Rx = simplejson.load(rx)
    with open('Ry_output.json', 'r') as ry:
        Ry = simplejson.load(ry)
    with open('Rz_output.json', 'r') as rz:
        Rz = simplejson.load(rz)
    with open('P_output.json', 'r') as p:
        P = simplejson.load(p)
    with open('Vx_output.json', 'r') as vx:
        Vx = simplejson.load(vx)
    with open('Vy_output.json', 'r') as vy:
        Vy = simplejson.load(vy)
    with open('Vz_output.json', 'r') as vz:
        Vz = simplejson.load(vz)  
    
    momentum_array = []
    for step in range(steps):
        system_momentum = 0
        for i in range(N):
            system_momentum += Vx[i][step]
            system_momentum += Vy[i][step]
            system_momentum += Vz[i][step]
        momentum_array.append(system_momentum)
    return momentum_array

def ang_momentum(N,Rx,Ry,Rz,Vx,Vy,Vz,m):
    momentx=[]
    momenty=[]
    momentz=[]
    ri=np.sqrt(Rx[0]**2+Ry[0]**2+Rz[0]**2)
    vi=np.sqrt(Vx[0]**2+Vy[0]**2+Vz[0]**2)
    momxi=m[i]*Vx[0]
    momyi=m[i]*Vy[0]
    momzi=m[i]*Vz[0]

    for i in xrange(N):
        r=np.sqrt(Rx[i]**2+Ry[i]**2+Rz[i]**2)
        v=np.sqrt(Vx[i]**2+Vy[i]**2+Vz[i]**2)
        theta=np.arccos(v/r)
        momx=m*Vx[i]*r*np.sin(theta)
        momy=m*Vy[i]*r*np.sin(theta)
        momz=m*Vz[i]*r*np.sin(theta)     
        
        momentx.append(momx-momxi)
        momenty.append(momy-momyi)
        momentz.append(momz-momzi)
        
    return momentx,momenty,momentz