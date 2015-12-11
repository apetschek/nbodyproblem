import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation
import simplejson
import numpy as np
from helper_functions import accel, masses

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

M_sun = 1.98 * 10**30 
M_param = M_sun 
m = masses(N,M_param)
steps = len(Rx[0])   

# Calculates total energy of the system
def energy(m,steps):

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

# Calculates total momentum of the system
def momentum(m,steps): 
    
    momentum_array = []
    
    for step in range(0,steps):
        
        system_momentum = 0
        
        for i in range(N):
            system_momentum += Vx[i][step]
            system_momentum += Vy[i][step]
            system_momentum += Vz[i][step]
        momentum_array.append(system_momentum)
    
    return momentum_array

def plot_energy(m,steps):
    E = energy(m,steps)
    plt.plot(E)
    plt.xlim(0,len(Rx[0]))
    plt.xlabel("Steps")
    plt.ylabel("System Energy (Joules)")
    plt.title("System Energy")
    plt.show()

def plot_momentum(m,steps):
    P = momentum(m,steps)
    plt.plot(P)
    plt.xlim(0,len(Rx[0]))
    plt.xlabel("Steps")
    plt.ylabel("System Momentum")
    plt.title("System Momentum")
    plt.show()

def plot_trajectory(Rx, Ry, N, limit):
    Rx = np.array(Rx) / AU
    Ry = np.array(Ry) / AU

    fig = plt.figure(figsize=(10,10))
    ax = fig.add_subplot(111, projection = '3d')
    
    for i in range(N):
        ax.plot(Rx[i], Ry[i], Rz[i])
    
    # ax.set_xlim(-limit, limit)
    # ax.set_ylim(-limit, limit)
    # ax.set_zlim(-limit, limit)
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    plt.show()

AU    = 1.496 * 10**10
limit = 1  #in AU

plot_trajectory(Rx, Ry, N, limit)
plot_energy(m,steps)
plot_momentum(m,steps)
