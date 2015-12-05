import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation
import simplejson
import numpy as np

with open('N_output.json', 'r') as n:
    N = simplejson.load(n)
with open('norms_output.json', 'r') as n:
    norms = simplejson.load(n)
with open('Rx_output.json', 'r') as rx:
    Rx = simplejson.load(rx)
with open('Ry_output.json', 'r') as ry:
    Ry = simplejson.load(ry)
with open('Rz_output.json', 'r') as rz:
    Rz = simplejson.load(rz)

def plot_trajectory(Rx, Ry, Rz, N, limit):
	Rx = np.array(Rx) / AU
	Ry = np.array(Ry) / AU
	Rz = np.array(Rz) / AU

	fig = plt.figure(figsize=(10,10))
	ax = fig.add_subplot(111, projection = '3d')
	
	for i in range(N):
		ax.plot(Rx[i], Ry[i], Rz[i])
	
	ax.set_xlim(-limit, limit)
	ax.set_ylim(-limit, limit)
	ax.set_zlim(-limit, limit)
	ax.set_xlabel('X')
	ax.set_ylabel('Y')
	ax.set_zlabel('Z')
	plt.show()

def plot_norms_hist(norms, N):
    fig, ((ax1, ax2),( ax3, ax4)) = plt.subplots(nrows=2, ncols=2)
    
    norms0 = np.array(norms[0]) / AU
    norms1 = np.array(norms[1]) / AU
    norms2 = np.array(norms[2]) / AU
    norms3 = np.array(norms[3]) / AU

    x_max = 2
    
    ax1.set_xlim(0, x_max)
    ax2.set_xlim(0, x_max)
    ax3.set_xlim(0, x_max)
    ax4.set_xlim(0, x_max)
    
    ax1.set_ylim(0, N)
    ax2.set_ylim(0, N)
    ax3.set_ylim(0, N)
    ax4.set_ylim(0, N)
    
    bins = [(x_max*i) / float(N) for i in range(N)]
    
    ax1.hist(norms0, bins=bins)
    ax2.hist(norms1, bins=bins)
    ax3.hist(norms2, bins=bins)
    ax4.hist(norms3, bins=bins)
    
    plt.show()

AU    = 1.496 * 10**10
limit = .1

plot_trajectory(Rx, Ry, Rz, N, limit)
plot_norms_hist(norms, N)