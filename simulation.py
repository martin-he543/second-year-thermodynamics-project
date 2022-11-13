"""
SIMULATION MODULE       Martin A. He
----------------------------------------------------------
A module for the simulation of the ball.
->
->

"""
import numpy as np
import ball as bl
import matplotlib.pyplot as plt

class Simulation:
    """
    SIMULATION CLASS
    Simulates the movement of hard spherical gas particles in a circular container.

    ATTRIBUTES
    - N_balls: the number of balls to simulate
    - r_balls: the radius of balls in the simulation
    - r_container: the radius of the container
    
    """
    def __init__(
        self,
        N_balls=1,
        r_balls=1,
        r_container=10,
    ):
        self.__N_balls=N_balls
        self.__r_balls=r_balls
        self.__r_container=r_container
        return

    def next_collision(self):
        """
        Sets up and performs the next collision.

        The logic of the next_collision method should be as follows:
        - find the time to the next collision
        - move the system to that point in time
        - perform the collision
        """

        return