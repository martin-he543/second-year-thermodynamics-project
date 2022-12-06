import numpy as np
import simulation as sim
import scipy.constants as spc
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import itertools as it
import heapdict as hd
import time as tm
import sys, os

#%% Simulation Plots ✔️
### SIMULATION PLOTTING METHODS
def plot_simulation_animation(m_balls_animation=5e-26, r_balls_animation=0.2,
                              N_balls_animation=100, collisions_animation=500,
                              r_container_animation=10,
                              random_speed_range_animation=500):
    sim_test_animation = sim.Simulation(m_balls=m_balls_animation, 
        r_balls=r_balls_animation, N_balls=N_balls_animation, r_container=
        r_container_animation, random_speed_range=random_speed_range_animation)

    sim_test_animation.run(collisions=collisions_animation, animate=True,
                           sim_title="1. Simulation")

def plot_distance(self, context="paper", absolute=True, relative=True, 
                    x_label_abs="Distance Distribution, Measured from O", 
                    y_label_abs= "Probability Density $/m^{-1}$",
                    plot_title_abs="Absolute Distance between Balls",
                    x_label_rel="Ball Distance, $/m$",
                    y_label_rel="Probability Density $/m^{-1}$",
                    plot_title_rel="Relative Distance, Measured between Balls",
                    style="darkgrid", palette="muted", font="sans-serif",
                    simulation_params=[5e-26,0.2,50,5000,10,500]):
    """ plot_distance | Plots the distance graphs for the simulation.
            PARAMETERS:
            -> context(str, optional):
            -> absolute(boolean, optional):
            -> relative(boolean, optional):
            -> x_label_abs(str, optional):
            -> y_label_abs(str, optional):
            -> plot_title_abs(str, optional):
            -> x_label_rel(str, optional):
            -> y_label_rel(str, optional):
            -> plot_title_rel(str, optional):
            -> style(str, optional):
            -> palette(str, optional):
            -> font(str, optional):
            -> simulation_params(list, optional):
    """
    _context, _style, _palette = context, style, palette
    _absolute, _relative = absolute, relative
    _x_label_abs, _y_label_abs = x_label_abs, y_label_abs
    _plot_title_abs = plot_title_abs
    _x_label_rel, _y_label_rel = x_label_rel, y_label_rel
    _plot_title_rel = plot_title_rel
    _font, _simulation_params = font, simulation_params

    m_balls_distance = _simulation_params[0]
    r_balls_distance = _simulation_params[1]
    N_balls_distance = _simulation_params[2]
    collisions_distance = _simulation_params[3]
    r_container_distance = _simulation_params[4]
    random_speed_range_distance = _simulation_params[5]

    sim_test = sim.Simulation(N_balls=N_balls_distance,
                            r_container=r_container_distance,
                            r_balls=r_balls_distance, m_balls=m_balls_distance,
                            random_speed_range=random_speed_range_distance)

    run_distance = sim_test.run(collisions=collisions_distance, KE=True,
                                distance_absolute=True, distance_relative=True)

    if _absolute: distance_absolute_run = run_distance["Absolute Distance"]
    if _relative: distance_relative_run = run_distance["Relative Distance"]

    if _absolute:
        print("Plotting Absolute Distance Distribution Graph...")
        plt.figure(num=_plot_title_abs)
        sns.set_theme(context=_context, style=_style, palette=_palette, font=_font)
        sns.distplot(distance_absolute_run)

        plt.title(_plot_title_abs)
        plt.xlabel(_x_label_abs); plt.ylabel(_y_label_abs)
        plt.tight_layout(); plt.show()

    if _relative:
        print("Plotting Relative Distance Distribution Graph...")
        plt.figure(num=_plot_title_rel)
        sns.set_theme(context=_context, style=_style, palette=_palette, font=_font)
        sns.distplot(distance_relative_run)

        plt.title(_plot_title_rel)
        plt.xlabel(_x_label_rel); plt.ylabel(_y_label_rel)
        plt.tight_layout(); plt.show()

def generate_dataset():
    """ generate_dataset | Generates a dataset
    """
    print("Generated")

def plot_ideal_gas_law():
    """ plot_ideal_gas_law | Plots the ideal gas law curves.
    """
    print("Generated")

def plot_van_der_waals():
    """ plot_van_der_waals | Plots the Van der Waals law curves.
    """
    print("Generated")

def plot_maxwell_boltzmann():
    """ plot_maxwell_boltzmann | Plots the Maxwell-Boltzmann distribution.
    """
    print("Generated")