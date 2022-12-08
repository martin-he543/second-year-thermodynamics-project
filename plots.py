import simulation as sim
import seaborn as sns
import matplotlib.pyplot as plt
import scipy as sp
import scipy.constants as spc
import concurrent.futures as ft
import numpy as np
import pandas as pd
import os, time

titleFont =     {'fontname': 'Kinnari', 'size': 13}
axesFont =      {'fontname': 'Kinnari', 'size': 9}
ticksFont =     {'fontname': 'SF Mono', 'size': 7}
errorStyle =    {'mew': 1, 'ms': 3, 'capsize': 3, 'color': 'blue', 'ls': ''}
pointStyle =    {'mew': 1, 'ms': 3, 'color': 'blue'}
lineStyle =     {'linewidth': 0.5}
lineStyleBold = {'linewidth': 1}
histStyle =     {'facecolor': 'green', 'alpha': 0.5, 'edgecolor': 'black'}

class Plots:
    def plot_simulation_animation(m_balls_animation=5e-26, r_balls_animation=0.2,
                                N_balls_animation=100, collisions_animation=500,
                                r_container_animation=10,
                                random_speed_range_animation=500):
        """ plot_simulation_animation | Plot the animation of the simulation.
                < PARAMETERS >
                -> m_balls_animation(str, optional):
                -> r_balls_animation

        """
        print("Commencing the Simulation Animation...")
        sim_test_animation = sim.Simulation(m_balls=m_balls_animation, 
            r_balls=r_balls_animation, N_balls=N_balls_animation, r_container=
            r_container_animation, random_speed_range=random_speed_range_animation)

        sim_test_animation.run(collisions=collisions_animation, animate=True)

    def generate_dataset(folder_name="data"):
        """ generate_dataset | Generates a dataset
            < PARAMETERS >
            -> generate
        """
        DATA_PATH = os.path.join(os.getcwd(), folder_name)
        if not os.path.exists(DATA_PATH): os.makedirs(DATA_PATH)
        m_ball = 5e-26
        l_N_ball = [10, 20]
        l_r_ball = [2, 4]
        l_r_container = [100, 200]
        l_random_speed_range = [500, 1000]
        collisions = 50

        print("Starting Simulations")

        for r_ball in l_r_ball:
            for N_ball in l_N_ball:
                for r_container in l_r_container:
                    for random_speed_range in l_random_speed_range:
                        fname = f"dataset_{N_ball}_{r_ball}_{r_container}_\
                                 {m_ball}_{random_speed_range}_{collisions}.csv"
                        FILE_PATH = os.path.join(DATA_PATH, fname)
                        if os.path.exists(FILE_PATH):
                            print(f"exists: {N_ball} balls, r_ball = {r_ball},\
                                    max speed = {random_speed_range},\
                                    r_container = {r_container},\
                                    {collisions} collisions")
                            continue
                        else:
                            s = sim.Simulation(N_balls=N_ball, r_container=r_container,
                                               r_balls=r_ball, m_balls=m_ball,
                                               random_speed_range=random_speed_range)
                            dataset = s.run(collisions=collisions, dataset=True)["dataset"]
                            dataset.to_csv(FILE_PATH)
                            print(f"Generated {FILE_PATH}")

Plots.plot_simulation_animation()
Plots.generate_dataset()