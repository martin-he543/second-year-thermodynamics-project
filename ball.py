"""
BALL MODULE       Martin A. He
----------------------------------------------------------
A module for the simulation of a bouncing ball. 
->
->

"""

import numpy as np
from copy import deepcopy

class Ball:
    def __init__(
        self,
        m_ball=1,
        r_ball=1,
        pos_ball=np.array(),
        vel_ball=np.array(),
    ):
        return

    def pos(self):
        """
        Return the current position.
        RETURNS 
            (np.ndarray w/ float values): the current position of the centre of the ball.
        """

    def vel(self):
        """
        Return the current velocity.
        RETURNS
            (np.ndarray w/ float values): the current 2D velocity of the ball.
        """

    def move(self,dt):
        """
        Move the ball to a new position: r' = r + v * dt.
        RETURNS

        """

    def time_to_collision(self,other):
        """
        Calculate the time until the next collision between this ball and another one, or the container.
        RETURNS
            
        """

    def collide(self,other):
        """
        Make the changes to the velocities of the ball, and the other one due to a collision.
        RETURNS

        """


