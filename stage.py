import pylab as pl
import matplotlib.pyplot as plt

# import simulation as sim
#%% Adjust your font settings
titleFont = {'fontname':'SF Mono','size':13}
axesFont = {'fontname':'SF Mono','size':9}
ticksFont = {'fontname':'SF Mono','size':7}
errorStyle = {'mew':1,'ms':3,'capsize':3,'color':'blue','ls':''}
lineStyle = {'linewidth':1}
radiusContainer = 10
xMin, xMax, yMin, yMax = -radiusContainer, radiusContainer, -radiusContainer, radiusContainer

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