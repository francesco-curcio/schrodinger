#!/usr/bin/env python3

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import animation
import subprocess
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('n_points', type=int)
parser.add_argument('n_t', type=int)
args = parser.parse_args()

subprocess.run("./main {} {}".format(args.n_points, args.n_t), shell=True)
t = np.genfromtxt("res.dat", delimiter=',', usecols=(0))
psi2 = np.genfromtxt("res.dat", delimiter=",")
psi2 = psi2[:, 1:-1]

max_psi2 = np.amax(psi2)

fig = plt.figure()
ax = plt.axes(xlim=(0, 1), ylim=(0, max_psi2))
line, = ax.plot([], [], lw=2)

# initialization function: plot the background of each frame
def init():
    line.set_data([], [])
    return line,

# animation function.  This is called sequentially
def animate(i):
    x = np.linspace(0, 1, args.n_points)
    y = psi2[i,:]
    line.set_data(x, y)
    return line,

# call the animator.  blit=True means only re-draw the parts that have changed.
anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=args.n_t, interval=20, blit=True)

# save the animation as an mp4.  This requires ffmpeg or mencoder to be
# installed.  The extra_args ensure that the x264 codec is used, so that
# the video can be embedded in html5.  You may need to adjust this for
# your system: for more information, see
# http://matplotlib.sourceforge.net/api/animation_api.html

anim.save('../../Desktop/particle.mp4', fps=60, extra_args=['-vcodec', 'libx264'])

plt.show()

