import numpy as np
import ball as bl
import matplotlib.pyplot as plt

class Event(tuple):
    """
    A tuple of 5 elements (ball_A, ball_B, count_A, count_B, dt).

    PARAMETERS
        ball_A (int): The first ball in impending collision.
        ball_B (int): The second ball in impending collision.
        count_A (int): The number of collisions the first ball 
            encountered prior to this impending collision calculation.
        count_B (int): The number of collisions the second ball 
            encountered prior to this impending collision calculation.
        dt (float): The global time this collision will happen on.
    """

    def ball_A(self):
        """
        RETURNS
            (int): The index of first ball in impending collision.
        """
        return self[0]

    def ball_B(self):
        """
        RETURNS
            (int): The index of second ball in impending collision.
        """
        return self[1]

    def count_A(self):
        """
        RETURNS
            (int): The number of collisions the first ball 
                encountered prior to this impending collision calculation.
        """
        return self[2]

    def count_B(self):
        """
        RETURNS
            (int): The number of collisions the second ball 
                encountered prior to this impending collision calculation.
        """
        return self[3]

    def dt(self):
        """
        RETURNS
            (float): The global time this collision will happen on.
        """
        return self[4]

    def pair(self):
        """
        RETURNS
            (list of int): A list of the two balls/ ball and container
                involved in collision.
        """
        return [self[0], self[1]]


class Simulation:
    """
    SIMULATION CLASS
    Simulates the movement of hard spherical gas particles in a circular container.

    PARAMETERS
    - N_balls: the number of balls to simulate
    - r_balls: the radius of balls in the simulation
    - r_container: the radius of the container
    
    """
    def __init__(
        self,
        N_balls=1,
        r_balls=1,
        m_balls=1,
        r_container=10
    ):
        self._N_balls=N_balls
        self._r_balls=r_balls
        self._r_container=r_container
        self._ball=[]
        self._temperature=[]
        self._KE=[]
        self._speed=[]
        self._brownian=[]
        
        self._distance_centre=[]

        for _ in range(0, N_balls):
            self._ball.append(bl.Ball(radius=r_balls, mass=m_balls))


    ### SIMULATION INFORMATION METHODS
    
    def __repr__(self):
        return ("Simulation properties: N = {self._N_ball} balls, r_ball = {self._r_ball}, m_ball = {self._m_ball}, r_container = {self._r_container}")
    def __str__(self):
        return ("Simulation: N = {self._N_ball} balls, r_ball = {self._r_ball, m_ball = {self._m_ball}, r_container = {self._r_container}")
    
    ### SIMULATION MOVEMENT METHODS
    
    def next_collision(self,ball,time):
        """
        Sets up and performs the next collision.

        The logic of the next_collision method should be as follows:
        - find the time to the next collision
        - move the system to that point in time
        - perform the collision
        """

    ### SIMULATION ATTRIBUTE METHODS
    
    def N_ball(self):
        """
        RETURNS
            (int): Number of balls in container.
        """
        return self._N_ball
    
    def ball(self):
        """
        RETURNS
            (list): A list of all the active balls in the simulation.
        """
        return self._ball
    
    def container(self):
        return self._container
    
    ### SIMULATION PROPERTY METHODS
    
    def pressure(self):
        return self._pressure

    def pressure_moyen(self):
        return self._pressure_moyen
    
    def temperature(self):
        return self._temperature
    
    def temperature_moyen(self):
        return self._temperature_moyen
    
    def KE(self):
        return self._KE
    
    def distance_centre(self):
        return self._distance_centre
    
    def speed(self):
        return self._speed
    
    def brownian(self):
        return self._brownian
    
    def dataset(self):
        return self._dataset
    
    ### SIMULATION OTHER METHODS
    
    def run(
        self,
        collisions = 10,
        time = 0.001
    ):
        """
        Runs the 2D simulation of colliding particles within the container.
        
        Arguments:
                collisions (int, optional): Number of collisions in the simulation.
                
        """