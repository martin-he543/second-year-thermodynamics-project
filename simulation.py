import numpy as np
import ball as bl
import matplotlib.pyplot as plt
import pandas as pd
import itertools as it
import heapdict as hd
import time as tm
import random

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
        N_balls = 1,                        # The number of balls.
        r_balls = 1,                        # The radius of the balls.
        m_balls = 1,                        # The mass of the balls.
        r_container = 10,                   # The radius of the container.
        
        random_position = True,             # Randomly-positioned balls.
        random_speed_range = 5              # Range for random speed generation.
    ):
        self._N_balls = N_balls             # The number of container balls.
        self._r_balls = r_balls             # The radius of the container balls.
        self._r_container = r_container     # The radius of the container.
        
        self._ball = []                     # All the balls of the container.
        self._temperature = []              # System temperature for all collisions.
        self._KE = []                       # System K.E. for all collisions.
        self._speed = []                    # Speed of all balls for all collisions.
        
        self._distance_centre = []          # Ball distance from origin.
        self._distance_relative = []        # Relative distance between balls.
        self._pairs = self.pair_combi()     # Lists all pair combinations of balls.
        self._random_position = random_position     # Sets random positioning.
        self._random_speed = random_speed_range     # Sets random speed range.
        
        self._brownian = []                 # Brownian motion investigation data.

        for _ in range(0, N_balls):
            self._ball.append(bl.Ball(radius=r_balls, mass=m_balls))
        if random_position == True:
            print("Random positioning enabled.")

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
        - Find the time to the next collision
        - Move the system to that point in time
        - Perform the collision
        """

    ### SIMULATION ATTRIBUTE METHODS
    
    def N_ball(self):
        """
        Returns how many container balls.
        RETURNS
            (int): Number of balls in container.
        """
        return self._N_ball
    
    def ball(self):
        """
        Lists all active balls in the simulation.
        RETURNS
            (list: bl.Ball): a list of all the active balls in the simulation.
        """
        return self._ball
    
    def container(self):
        """
        Provides the container of the simulation.
        RETURNS
            (bl.Container): the simulation's container object.
        """
        return self._container
    
    ### SIMULATION PROPERTY METHODS
    
    def pressure(self):
        """
        Provides the pressure values for every __ container collisions.
        N.B: Must enable Simulation.run(pressure = True).
        RETURNS
            pandas shit
            pressure(float): system pressure.
        """
        return self._pressure

    def pressure_moyen(self):
        """
        Provides the average pressure of the system.
        N.B: Must enable Simulation.run(pressure = True).
        RETURNS
            pressure_moyen(float): the system's average steady-state pressure.
        """
        return self._pressure_moyen
    
    def temperature(self):
        """
        Provides the system temperature at all collision times.
        N.B: Must enable Simulation.run(temperature = True).
        RETURNS
        """
        return self._temperature
    
    def temperature_moyen(self):
        """
        Provides the average temperature of the system.
        N.B: Must enable Simulation.run(temperature = True).
        RETURNS
            temperature_moyen(float): the system's average temperature.
        """
        return self._temperature_moyen
    
    def KE(self):
        """
        Provides the total system KE for all collisions.
        N.B: Must enable Simulation.run(KE = True).
        RETURNS
        """
        return self._KE
    
    def distance_centre(self):
        """
        Provides distances from the origin for all balls, in all collisions.
        N.B: Must enable Simulation.run(distance_centre = True).
        RETURNS
            distance_centre(list of (float)): ball distances from the origin for all collisions.
        """
        return self._distance_centre
    
    def distance_relative(self):
        """
        Provides relative distance between all balls, in all collisions.
        N.B: Must enable Simulation.run(distance_relative = True).
        RETURNS
            relative(list of (float)): ball distances between all possible pairs of balls, in all collisions.
        """
        return self._distance_centre
    
    def speed(self):
        """
        Provides the speed for all balls, in all collisions.
        N.B: Must enable Simulation.run(speed = True).
        RETURNS
            speed(list of (float)): ball speeds, in all collisions.
        """
        return self._speed
    
    def brownian(self):
        """
        Provides the Brownian motion investigation dataset.
        N.B: Must enable Simulation.run(brownian = True).
        RETURNS
            pandas shit
        """
        return self._brownian
    
    def dataset(self):
        """
        Provides a complete dataset, for this simulation.
        N.B: Must enable Simulation.run(dataset = True).
        RETURNS
        """
        return self._dataset
    
    ### SIMULATION RANDOMISATION METHODS
    
    def gen_random_positions(self, start=0):
        """
        Generates the non-overlapping random ball positions.
        
        PARAMETERS
            start (boolean, optional):
        RAISES
            Exception:
                Balls cannot fit in this container. Please reduce N_balls or increase r_container to proceed.
        """
        def gen_random_uniform(maximum_range):
            """
            Provides a uniform distribution centered at (0,0), generating random floats.
            PARAMETERS
                maximum_range (float): sets the maximum range for the distribution.
            RETURNS
                (float): a random float between [-max_range, max_range].
            """
            return random.uniform(-maximum_range,maximum_range)
        
        for i in range(start, self._N_balls):
            position, error_count = np.zeros(2), 0
            
            while True:
                # Extreme case error handling to prevent computational overload.
                if error_count > 1e5:
                    raise Exception \
                        ("The area of this container is too small for ball size.")
                x = gen_random_uniform(self._r_container - self._ball[i]._radius)
                y = gen_random_uniform(self._r_container - self._ball[i].radius)
                # Check if the randomly-assigned position is valid.
                while(np.sqrt(x**2 + y**2) >= self._r_container - self._ball[i]._radius):
                    x = gen_random_uniform(self._r_container - self._ball[i]._radius)
                    y = gen_random_uniform(self._r_container - self._ball[i].radius)
                position = np.array([x,y])
                append = False
                
                for j in range(0, i):
                    distance = np.sqrt((self._ball[j]._pos_ball[0] - position[0])**2 \
                        + (self._ball[j]._pos_ball[1] - position[1]**2))
                    if distance <= self._ball[i].radius + self._ball[j]._radius:
                        append = False
                        error_count += 1
                        break
                    else: append = True
                if append or i == 0: break
            self.ball[i].set_pos(position)

    ### SIMULATION BROWNIAN MOTION INVESTIGATION
    
    def brownian_init(self, radius = 5, mass = 10):
        """
        Initialisation of Brownian motion investigation, ball 0.
        -> Position of ball 0 is initialised to (0,0).
        -> Remaining balls are randomly distributed in the container.
        PARAMETERS
            radius (float, optional): radius of ball being investigated.
            mass (float, optional): mass of ball being investigated.
        """
        self._ball[0].set_pos(np.array([0.0,0.0]))
        self._ball[0].set_radius(radius)
        self._ball[0].set_mass(mass)
        self.gen_random_positions(start = 1)

    def set_ball_vel(self, velocity_list):
        """
        Sets the velocity of all the balls from a list of velocities.
        PARAMETERS
            velocity_list (list of (np.ndarray of (floats))): lists all the ball velocities in their x- and y- directions.
        """
        for i, velocity in enumerate(velocity_list):
            self._ball[i].set_vel(velocity)

    def brownian_velocities(self, maximum_speed):
        """
        Generates and assigns random velocities to all balls from a uniform random distribution of x, and y-components of velocity.
        PARAMETERS
            maximum_speed (float): the range to generate velocities from.
        """
                
        def gen_random_velocities(N_balls, maximum_speed_range):
            """
            Generates random velocities from ma uniform distribution for a given number of balls in the range of [-maximum_speed_range, maximum_speed_range].
            PARAMETERS
                N_balls (int): the number of balls.
                random_speed_range(float): the ± range in the x, and y-components.
            RETURNS
            """
            def gen_random_uniform(maximum_range):
                """
                Provides a uniform distribution centered at (0,0), generating random floats.
                PARAMETERS
                    maximum_range (float): sets the maximum range for the distribution.
                RETURNS
                    (float): a random float between [-max_range, max_range].
                """
                return random.uniform(-maximum_range,maximum_range)
            
            list = []
            for _ in range(N_balls):
                list.append(np.array(gen_random_uniform(maximum_speed_range),\
                    gen_random_uniform(maximum_speed_range)))
            return list
                    
        list = gen_random_velocities(self._N_balls, self._random_speed)
        self.set_ball_vel(list)

        
    ### SIMULATION OTHER METHODS
    
    ### SIMULATION RUN METHODS
    
    def run(
        self,
        collisions = 10,            # The number of collisions.
        time = 0.001,               # Time period between collisions.
        
        pressure = False,           # Enable pressure calculations.
        temperature = True,         # Enable temperature calculations.
        KE = False,                 # Enable kinetic energy calculations.
        speed = False,              # Enable speed calculations.
        brownian = False,           # Enable Brownian motion inveestigation.
        dataset = False,            # Enable complete dataset.
        test_pressure = False,
        test_temperature = False,
        distance_centre = False,    # Use ball distance from the origin.
        distance_relative = False,  # Use relative distances between balls.
        progress_bar = True         # Enable progress bar animation in terminal.
    ):
        """
        Runs the 2D simulation of colliding particles within the container.
        
        PARAMETERS
        -> collisions (int, optional): number of collisions in the simulation.
        -> time (float, optional): time period (s) between animation frames.
            
        -> pressure (boolean, optional): records pressure for every __ collisions with the wall of the container.
        -> temperature (boolean, optional): records temperature.
        -> KE (boolean, optional): records the system KE for every collision.
        -> speed (boolean, optional): records speed of all balls in all collisions.
        -> brownian (boolean, optional): records data for Brownian motion.
        -> dataset (boolean, optional): records dataset of simulation information.
        -> distance_centre (boolean, optional): records the distances of all balls from the origin, in all collisions.
        -> distance_relative (boolean, optional): records the relative distances between all the balls, in all collisions.
        -> progress_bar (boolean, optional): displays a progress bar.
            
        RETURNS
        
        """
        

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
            (int): the index of first ball in impending collision.
        """
        return self[0]

    def ball_B(self):
        """
        RETURNS
            (int): the index of second ball in impending collision.
        """
        return self[1]

    def count_A(self):
        """
        RETURNS
            (int): the number of collisions the first ball encountered prior to this impending collision calculation.
        """
        return self[2]

    def count_B(self):
        """
        RETURNS
            (int): the number of collisions the second ball encountered prior to this impending collision calculation.
        """
        return self[3]

    def dt(self):
        """
        RETURNS
            (float): the global time this collision will happen on.
        """
        return self[4]

    def pair(self):
        """
        RETURNS
            (list of int): a list of the two balls / ball and container involved in collision.
        """
        return [self[0], self[1]]