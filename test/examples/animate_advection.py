import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

field  = np.load('Advection.npy')
coords = np.load('coordinates.npy')
times  =np.load('times.npy')

fig, ax = plt.subplots(figsize=(5, 3))
ax.set(xlim=(0, 9), ylim=(-2, 2))

def animate(i):
    ax.set_title('Frame %i' %i)
    plt.plot(coords[0,:], field[i,0,:])
    plt.plot(coords[1,:], field[i,1,:])
    plt.plot(coords[2,:], field[i,2,:])

anim = FuncAnimation(
    fig, 
    animate, 
    interval=300, 
    frames=len(times)-1,
    repeat=False
    )

plt.draw()
plt.show()
