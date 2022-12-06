import simulation as sim
import numpy as np
import scipy.constants as spc
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import itertools as it
import heapdict as hd
import time as tm

m_balls_animation = 5e-26
r_balls_animation = 0.2
N_balls_animation = 100
collisions_animation = 500
r_container_animation = 10
random_speed_range_animation = 500


sim_test_animation = sim.Simulation(N_balls=N_balls_animation, r_container=\
    r_container_animation, r_balls=r_balls_animation, m_balls=m_balls_animation,
    random_speed_range=random_speed_range_animation)
run = sim_test_animation.run(collisions=collisions_animation, animate=True,
                                sim_title="1. Simulation")

print("Plotting Graph 3 of 4")

plt.figure(num="Relative Distance between Balls")
sns.set(context="paper", style="darkgrid", palette="muted")

sns.displot(dist_rel_test)

plt.title("Relative Ball Distance Distribution")
plt.xlabel(r"Distance between Balls $/m$")
plt.ylabel(r"Probability Density $/m^{-1}$")
plt.tight_layout()
plt.show()


print("Plotting Graph 4 of 4")

plt.figure(num="Distance of Balls from Origin")
sns.set(context="paper", style="darkgrid", palette="muted")

sns.distplot(dist_centre_test, kde=False, norm_hist=True)

plt.title("Distance Distribution from Origin")
plt.xlabel("Ball Distance from Origin $/m$")
plt.ylabel("Probability Density $/m^{-1}$")
plt.tight_layout()
plt.show()

print("End of Script")