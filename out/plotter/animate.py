"""
A simple example of an animated plot
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

fig, ax = plt.subplots()

plt.ylim([-.6,1.2])

x = np.arange(-20., 20, 0.01)
#line, = ax.plot(x, np.sin(x))

data = np.genfromtxt("out/data/Evo.dat", delimiter=',')

xt = data[:100,1]
tt = data[::100,0]

#X, Y = np.meshgrid(et,xt)
ENE = data[:,2].reshape(len(tt),len(xt))
VEL = data[:,3].reshape(len(tt),len(xt))

print(x)

line, = ax.plot(xt, ENE[0])
ax.plot(xt, ENE[0])

def animate(i):
    line.set_ydata(ENE[i])
    return line,


# Init only required for blitting to give a clean slate.
def init():
    line.set_ydata(np.ma.array(xt, mask=True))
    return line,

ani = animation.FuncAnimation(fig, animate, np.arange(1, 79), init_func=init,
                              interval=20, blit=True)
plt.show()
