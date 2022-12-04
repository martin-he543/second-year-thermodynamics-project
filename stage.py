import simulation as sm
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
import pandas as pd
import concurrent.futures as ft
import os, time


#%%  
# Task 1 
# GENERATE CSV DATASETS for later reference

DATA_PATH = os.path.join(os.getcwd(), "data")
if not os.path.exists(DATA_PATH): os.makedirs(DATA_PATH)

m_ball = 5e-26
l_N_ball = [10, 20]; l_r_ball = [2,4]; l_r_container = [100,200]
l_random_speed_range = [500,1000]; collisions = 50

print("Simulations beginning.")
for r_ball in l_r_ball:
    for N_ball in l_N_ball:
        for r_container in l_r_container:
            for random_speed_range in l_random_speed_range:
                fname = f"dataset_{N_ball}_{r_ball}_{r_container}_{m_ball}\
                    _{random_speed_range}_{collisions}.csv"
                FILE_PATH = os.path.join(DATA_PATH,fname)
                
                if os.path.exists(FILE_PATH):
                    print(f"exists: {N_ball} balls, r_ball = {r_ball}, \
                            max speed = {random_speed_range}, \
                            r_container = {r_container}, {collisions} collisions")
                    continue
                else:
                    s = sm.Simulation(
                        N_balls=N_ball,
                        r_container=r_container,
                        r_balls=r_ball,
                        m_balls=m_ball,
                        random_speed_range=random_speed_range)
                    dataset = s.run(collisions=collisions, dataset=True)["dataset"]
                    dataset.to_csv(FILE_PATH)
                    print(f"Generated {FILE_PATH}")
print("End.")