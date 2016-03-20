"""
A simple example of an animated plot
"""
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation

fig, ax = plt.subplots()

plt.ylim([-.6,1.2])

x = np.arange(0, 40, 0.01)
#line, = ax.plot(x, np.sin(x))

data = np.genfromtxt("out/data/Evo.dat", delimiter=',')

xt = data[:100,1]
yt = data[::100,0]

X, Y = np.meshgrid(xt,yt)
Z = data[:,2].reshape(len(yt),len(xt))

line, = ax.plot(xt, Z[0])

#print(xt)


def animate(i):
    line.set_ydata(Z[i])
    return line,


# Init only required for blitting to give a clean slate.
def init():
    line.set_ydata(np.ma.array(xt, mask=True))
    return line,

ani = animation.FuncAnimation(fig, animate, np.arange(1, 200), init_func=init,
                              interval=20, blit=True)
plt.show()
