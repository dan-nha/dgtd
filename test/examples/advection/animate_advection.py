import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

field  = np.load('Advection.npy')
coords = np.load('coordinates.npy')
times  =np.load('times.npy')

fig, ax = plt.subplots(figsize=(5, 3))
ax.set(xlim=(0, 9), ylim=(-1.5, 1.5))

for i in range(len(times)-1):
    ax.clear()        

    plt.plot(coords[0,:], field[i,0,:])
    plt.plot(coords[1,:], field[i,1,:])
    plt.plot(coords[2,:], field[i,2,:])

    plt.pause(0.01)

plt.show()
plt.clf() 
