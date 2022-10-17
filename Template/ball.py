"""
BALL MODULE       Martin A. He
----------------------------------------------------------
A module for the simulation of 
->
->

"""

import numpy as np
from copy import deepcopy

class Ball:
    """
    BALL CLASS
    ________

    Attributes:
        mass (float): The mass of the ball.
        radius (float): The radius of the ball.
        pos (numpy.ndarray of float): The position of the ball.
        vel (numpy.ndarray of float): The velocity of the ball.
        count (int): The number of collisions the ball experienced.
    """
    
    def __init__(
        self,
        N_ball=1,
        r_container=10,
        r_ball=1,
        m_ball=1,
        random_pos=True,
        random_speed_range=5,
    ):
        return