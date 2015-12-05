import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.animation as animation
import simplejson

with open('test_N_output.json', 'r') as n:
    N = simplejson.load(n)
with open('test_Rx_output.json', 'r') as rx:
    Rx = simplejson.load(rx)
with open('test_Ry_output.json', 'r') as ry:
    Ry = simplejson.load(ry)
with open('test_Rz_output.json', 'r') as rz:
    Rz = simplejson.load(rz)

def plot_output(Rx, Ry, Rz, N, limit):

	fig = plt.figure(figsize=(10,10))
	ax = fig.add_subplot(111, projection = '3d')
	for i in range(N):
		ax.plot(Rx[i], Ry[i], Rz[i])
	ax.set_xlim(-limit,limit)
	ax.set_ylim(-limit,limit)
	ax.set_zlim(-limit,limit)
	ax.set_xlabel('X')
	ax.set_ylabel('Y')
	ax.set_zlabel('Z')
	plt.show()

AU    = 1.496 * 10**10
limit = .1*AU
plot_output(Rx, Ry, Rz, N, limit)