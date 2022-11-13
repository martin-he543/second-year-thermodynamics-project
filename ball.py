import numpy as np
from copy import deepcopy

class Ball:
    n_ball = 0
    def __init__(
        self,
        m_ball=1,
        r_ball=1,
        pos_ball=np.array([0,0]),
        vel_ball=np.array([0,0]),
    ):
        self.__m_ball = m_ball
        self.__r_ball = r_ball
        self.__pos_ball = pos_ball
        self.__vel_ball = vel_ball
        n_ball += 1
        return

    def pos(self):
        """
        Return the current position.
        RETURNS 
            (np.ndarray w/ float values): 
            the current position of the centre of the ball.
        """
        return self.__pos_ball

    def vel(self):
        """
        Return the current velocity.
        RETURNS
            (np.ndarray w/ float values): 
            the current 2D velocity of the ball.
        """
        return self.__vel_ball

    def move(self,dt):
        """
        Move the ball to a new position: 
        r' = r + v * dt.
        """
        self.__pos_ball = self.__pos_ball + (self.__vel_ball * dt)

    def time_to_collision(self,r,v,R):
        """
        Calculate the time until the next collision between this ball and another one, or the container.
        """
        return -(np.dot(r,v)/np.dot(v,v)) + np.sqrt((np.dot(r,v))**2 - np.dot(v,v)(np.dot(r,r) - R*22))

    def collide(self,other):
        """
        Make the changes to the velocities of the ball, and the other one due to a collision.
        RETURNS

        """
        

