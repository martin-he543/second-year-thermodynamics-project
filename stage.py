import simulation as sim
import matplotlib.pyplot as plt
import seaborn as sns

m_balls_animation = 5e-26
r_balls_animation = 0.2
N_balls_animation = 100
collisions_animation = 500
r_container_animation = 10
random_speed_range_animation = 500

m_balls_distance = 5e-26
r_balls_distance = 0.2
N_balls_distance = 50
collisions_distance = 5000
r_container_distance = 10
random_speed_range_distance = 500

#%% Stage Animation
### STAGE ANIMATION
sim_test_animation = sim.Simulation(m_balls=m_balls_animation, 
    r_balls=r_balls_animation, N_balls=N_balls_animation, r_container=
    r_container_animation, random_speed_range=random_speed_range_animation)

run_animation = sim_test_animation.run(collisions=collisions_animation, 
                                       animate=True, sim_title="1. Simulation")

#%% Stage Distance Calculations
### STAGE DISTANCE CALCULATIONS

sim_test = sim.Simulation(N_balls=N_balls_distance,
                          r_container=r_container_distance,
                          r_balls=r_balls_distance, m_balls=m_balls_distance,
                          random_speed_range=random_speed_range_distance)

run_distance = sim_test.run(collisions=collisions_distance, KE=True,
                            distance_absolute=True, distance_relative=True)

dist_centre_test = run_distance["distance from centre"]
dist_rel_test = run_distance["relative distance"]



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