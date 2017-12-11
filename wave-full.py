#--- EDP TP MSIAM 1
#--- TP Wave equation
#--- Maelle Nodet, 2016
# Acknowledgement for the animation: author: Jake Vanderplas
# email: vanderplas@astro.washington.edu website: http://jakevdp.github.com

import numpy as np
import scipy
import math
from matplotlib import pyplot as plt
from matplotlib import animation

#--- Initial set up

# independant parameters
n = 90 # grid point in space
dx = 1 / n # step in space
x = np.arange(0, 1+dx, dx)# [0:dx:1]
sx = x.size 
m = 100 # nb of time steps
T = 3
dt = 1 / m
t = np.arange(0, T+dt, dt)# 0:dt:T
st = t.size
c = 1

# initial conditions 
# f = lambda z: np.sin(math.pi * z)
f = lambda z: (np.logical_and((z > 0), (z < 1))) * np.sin(5. * math.pi * z)
g = lambda z: 0 * z
u0 = f(x)
v0 = g(x)
# plt.figure()
# plt.plot(x, f(x),label="f")
# plt.plot(x, g(x),label="g")
# plt.legend()
# plt.show()

#--- PDE numerical resolution : first question

# Define A, used for the first scheme
A = np.zeros((sx, sx))
for k in range(1, sx-1):
    A[k, k - 1] = (c*dt/dx)**2 
    A[k, k    ] = 2 * (1. - (c*dt/dx) ** 2)
    A[k, k + 1] = (c*dt/dx)**2
# boundary conditions
A[0,    0   ] = 1
A[sx-1, sx-1] = 1
# Define B
B = -1. * np.identity(sx) 

# Initialisation: u0
u = np.zeros((st, sx))
u[0,:] = u0

# First iteration: u1, order 1
u[1,:] = u[0,:] + dt * v0

# numerical 2-step scheme
for i in range(1, st-1):
	u[i+1,:]    = A.dot(u[i,:]) + B.dot(u[i-1,:])
	u[i+1,0]    = 0. # boundary condition
	u[i+1,sx-1] = 0. # boundary condition
	#plt.plot(x,u[i+1,:])
	#plt.show()

# animation for u
# First set up the figure, the axis, and the plot element we want to animate
fig = plt.figure(2)
ax = plt.axes(xlim=(0, 1), ylim=(-1, 1))
line, = ax.plot([], [])

# initialization function: plot the background of each frame
def init():
    line.set_data([], [])
    return line,

# animation function.  This is called sequentially
def animate(i):
    line.set_data(x, u[i,:])
    return line,

# call the animator.  blit=True means only re-draw the parts that have changed.
# careful, on a Mac blit=True is a known bug => set it to False!
anim = animation.FuncAnimation(fig, animate, init_func=init,
                               frames=st, interval=100, blit=False)
plt.show()

# COMPLETE FOR THE REST OF THE PRACTICAL SESSION QUESTIONS

# IF YOU CHOOSE THIS SESSION FOR YOUR FINAL EXAM:
# REACHABLE GRADES:
# If you stop at question 5: 
#		- acceptable: if report is reasonably good, and most questions are ok
#		- exceeds expectations: if report is great and all questions are correct
#
# If you stop at question 7:
#		- acceptable: if report is acceptable and some questions are ok
#		- exceeds expectations: if report is good and most questions are correct
#		- outstanding: great report, everything is correct
#
# If you also do the 2D question (q 8):
# 		- exceeds expectations: if report is acceptable and most questions 1-7 are correct, but question 8 does not work well
#		- outstanding: good report, most questions are correct, including question 8



