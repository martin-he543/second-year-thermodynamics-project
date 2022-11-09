import numpy as np
import scipy as sp
import pylab as pl
import matplotlib.pyplot as plt
import matplotlib.animation as anim
from copy import deepcopy

#%% Adjust your font settings
titleFont = {'fontname':'SF Mono','size':13}
axesFont = {'fontname':'SF Mono','size':9}
ticksFont = {'fontname':'SF Mono','size':7}
errorStyle = {'mew':1,'ms':3,'capsize':3,'color':'blue','ls':''}
lineStyle = {'linewidth'
             :1}

radiusContainer = 10
xMin, xMax, yMin, yMax = -radiusContainer, radiusContainer, -radiusContainer, radiusContainer

# #%% ========================= BALL CLASS FILE ==============================
# Class for the "ball.py" file.
# class Ball:
#     def _init_( 
#         # Initialises the existence of the ball.
#         self, m_ball=1, r_ball=1,
#         pos_ball = np.array(),
#         vel_ball = np.array(),
#     ):
#         self.__m_ball = m_ball
#         self.__r_ball = r_ball
#         self.__pos_ball = pos_ball
#         self.__vel_ball = vel_ball
#         return

#     def pos(self): 
#         # Returns the current position of the centre of the ball.
#         return self.__pos_ball

#     def vel(self): 
#         # Returns the current velocity of the ball.
#         return self.__vel_ball

#     def move(self,dt): 
#         # Moves the ball to a new position r' = r + v * dt
#         self._pos_ball = self.pos_ball + (self._vel_ball * dt)
#         return

#     def time_to_collision(self,other): 
#         # Calculates the time until the next collision between this ball and another one, or the container.
#         return
        
#     def collide(self,other): 
#         # Makes the changes to the velocities of the ball and the other due to a collision.
#         return

#%% ========================= STAGE CLASS FILE ==============================
# Class for the "stage.py" file.
firstPlot = pl.figure()

ax = pl.axes(xlim=(xMin, xMax), ylim=(yMin, yMax))
patch1 = pl.Circle([0., 0.], 1, ec='red', fc='red')
patch2 = pl.Circle([6., 6.], 1, ec='blue', fc='blue')
ax.add_patch(patch1)
ax.add_patch(patch2)
container = pl.Circle([0., 0.], 10, ec='black', fc='none')
ax.add_patch(container)

plt.title('Thermodynamics Snookered', **titleFont)
plt.xlabel("x Position / m",**axesFont)
plt.ylabel("y Position / m",**axesFont)
plt.xticks(**ticksFont)
plt.yticks(**ticksFont)
plt.gca().set_aspect('equal')

# How to Move Objects
for i in range(-10, 10):
    patch1.center = [i, i]
    pl.pause(0.001)
pl.show()


#%% ========================= SIMULATION CLASS FILE ==============================
# Class for the "simulation.py" file.


# class Simulation:
#     def _init_(
#         self, N_balls=1, r_balls=1,
#         r_container=10,
#     ):
#         self._N_balls=N_balls, self.r_balls=r_balls, self._r_container=r_container
#         return

#     def next_collision(self):
#         return