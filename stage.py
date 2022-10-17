import pylab as pl
import simulation as sim


f = pl.figure()
patch1 = pl.Circle([0., 0.], 4, fc='r')
patch2 = pl.Circle([5., 2.], 4, fc='b')
ax = pl.axes(xlim=(-10, 10), ylim=(-10, 10))
ax.add_patch(patch1)
ax.add_patch(patch2)
pl.show()