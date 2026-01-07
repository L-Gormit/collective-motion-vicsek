import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation

# Load data
data = np.loadtxt("Traj.txt")

# Columns
t  = data[:,0]
i  = data[:,1].astype(int)
x  = data[:,2]
y  = data[:,3]
th = data[:,4]

# Simulation parameters
N = i.max() + 1
L = 10.0

# Reshape by time
times = np.unique(t)
nt = len(times)

X = x.reshape(nt, N)
Y = y.reshape(nt, N)
TH = th.reshape(nt, N)

# Set up figure
fig, ax = plt.subplots()
ax.set_xlim(0, L)
ax.set_ylim(0, L)

scat = ax.scatter(X[0], Y[0], s=10)

# Velocity arrows
quiv = ax.quiver(X[0], Y[0],
                 np.cos(TH[0]), np.sin(TH[0]),
                 scale=50)

def update(frame):
    scat.set_offsets(np.c_[X[frame], Y[frame]])
    quiv.set_offsets(np.c_[X[frame], Y[frame]])
    quiv.set_UVC(np.cos(TH[frame]), np.sin(TH[frame]))
    ax.set_title(f"t = {times[frame]:.2f}")
    return scat, quiv

ani = FuncAnimation(fig, update, frames=nt, interval=50)

plt.show()

