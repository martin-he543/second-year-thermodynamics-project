    #%% Simulation Initialisation
"""
SIMULATION MODULE (import simulation as sm) | Created by 𝑀𝒶𝓇𝓉𝒾𝓃 𝒜. 𝐻𝑒, 2022.11.20
For the simulation of elastic collisions of balls with others, and container.
"""
import numpy as np
import ball as bl
import matplotlib.pyplot as plt
import pandas as pd
import itertools as it
import heapdict as hd
import time as tm
import scipy.constants as spc
import sys, os

class Simulation:
    """
        SIMULATION CLASS | 
        Simulates the movement of hard spherical gas particles in a circular container.
        PARAMETERS
            N_balls: the number of balls to simulate
            r_balls: the radius of balls in the simulation
            r_container: the radius of the container
    
    """
    def __init__(self,
        N_balls = 1,                        # The number of balls.
        r_balls = 1,                        # The radius of the balls.
        m_balls = 1,                        # The mass of the balls.
        r_container = 10,                   # The radius of the container.
        
        random_position = True,             # Randomly-positioned balls.
        random_speed_range = 5              # Range for random speed generation.
    ):
        self._N_balls = N_balls             # The number of container balls.
        self._m_balls = m_balls             # The mass of the container balls.
        self._r_balls = r_balls             # The radius of the container balls.
        self._r_container = r_container     # The radius of the container.
        
        self._temperature = []              # System T, for all collisions.
        self._KE = []                       # System K.E...
        self._speed = []                    # Speed of all balls...
        self._balls = []                    # List of all balls.
        
        self._distance_absolute = []        # Ball distance from origin.
        self._distance_relative = []        # Relative distance between balls.
        self._pairs = self.pair_combn()     # Lists all pair combinations of balls.
        
        self._random_position = random_position # Sets random positioning.
        self._random_speed = random_speed_range # Sets random speed range.
        self._brownian = []                 # Brownian motion investigation data.
        
        self._container = bl.Container(radius=r_container) # Container choice.
        self._N_container_collisions = 0    # Number of container collisions.
        self._dp_container = []             # Changes in container momentum.
        
        self._pq = hd.heapdict()            # Priority queueing.
        self._global_time = 0               # Global timer.
        self._events = []                   # List of current events.
        self._min_dt = 0                    # Minimum time to next collision.
        self._N_collisions = 0              # Number of collisions.

        for _ in range(0, N_balls):
            self._balls.append(bl.Ball(radius=r_balls, mass=m_balls))
        if random_position: self.gen_random_positions()    
        self.gen_random_velocities(maximum_speed_range=random_speed_range)

    #%% Simulation Information
    ### SIMULATION INFORMATION METHODS
    # Gives all the simulation methods for returning information on the Simulation.
    def __repr__(self):
        return ("Simulation properties: N = {self._N_balls} balls, r_balls = {self._r_balls}, m_balls = {self._m_balls}, r_container = {self._r_container}")
    
    def __str__(self):
        return ("Simulation: N = {self._N_balls} balls, r_balls = {self._r_balls}, m_balls = {self._m_balls}, r_container = {self._r_container}")
    
    def glossary(self, speed=False, KE=False,
                 distance_absolute=False, distance_relative = False,temperature=False, test_temperature=False,
                 pressure=False, test_pressure=False,
                 dataset=False, brownian=False):
        """ glossary | Appends glossary of various datasets at end of simulation.
        """
        gloss = {}
        if speed: gloss["Speed"] = self._speed
        if KE: gloss["Kinetic Energy"] = self._KE
        if distance_absolute: gloss["Distance from O"] = self._distance_absolute
        if distance_relative: gloss["Relative Distance"] = self._distance_relative
        if temperature: gloss["<Temperature>"] = self._temperature_moyen
        if test_temperature: gloss["Temperature"] = self._temperature
        if pressure: gloss["<Pressure>"] = self._pressure_moyen
        if test_pressure: gloss["Pressure"] = self._pressure
        if dataset: gloss["Dataset"] = self._dataset
        if brownian: gloss["Brownian"] = self.brownian
        return gloss
    
    #%% Simulation Movement
    ### SIMULATION MOVEMENT METHODS
    # Gives all the simulation methods for defining movement of balls.
    # For clarity, next_collision function has been split into several subfunctions.
    """ next_collision | Sets up and performs the next collision.
        The logic of the next_collision method should be as follows:
        1. Find the time to the next collision
        2. Move the system to that point in time
        3. Perform the collision. 
    """
    def init_collision_time(self):
        """ init_collision_time | Initialise first next collision time calculation.
            Calculate all possible ball pairs, and their impending collision times.
            Collision times are recorded as "Event" objects (see Event).
            All collision events are added.
        """
        for pair in self._pairs:    # All possible combos of ball pairs.
            ball_A, ball_B = self._balls[pair[0]], self._balls[pair[1]]
            dt = ball_A.time_to_collision(ball_B)
            if dt != np.inf:        # Only considers possible collisions.
                self._pq[Event((pair[0],pair[1],ball_A._count,ball_B._count,dt))] = dt
        
        for i, ball in np.ndenumerate(self._balls): # Change to enumerate if fails
            dt = ball.time_to_collision(self._container)
            if dt != np.inf:                    # Consider possible collisions.
                 self._pq[Event((i,self._N_balls,ball._count,-1,dt))] = dt

    def collision_time(self):
        """ collision_time | Finds next times of collision for balls that collided. """
        collided_ball = set()
        for event in self._events:
            for collided in event.pair():
                collided_ball.add(collided)
        
        for element in collided_ball:       # Add events to the priority queue.
            if element != self._N_balls:    # Collisions with the container.
                dt = self._balls[element].time_to_collision(self.container)
                if dt != np.inf:
                    self._pq[Event((element,self._N_balls,self._balls[element]._count,\
                        -1,dt+self._global_time,))] = (dt + self._global_time)
        
                for j in range(self._N_balls):  # Collisions with other balls.
                    if j != element:
                        if j < element:         # Ensure smaller index first.
                            ball_A, ball_B = self._balls[j],self._balls[element]
                            index_A, index_B = j, element
                        else:
                            ball_A, ball_B = self._balls[element],self._balls[j]
                            index_A,index_B = element, j
                        
                        dt = ball_A.time_to_collision(ball_B)
                        if dt != np.inf:        # Check collision is possible.
                            self._pq[Event((index_A,index_B,self._balls[index_A]._count,self._balls[index_B]._count,dt + self._global_time))]\
                                = (dt + self._global_time)
                                
            else:   # Container undergoes a collision.
                for j in range(self._N_balls):
                    dt = self._balls[j].time_to_collision(self._container)
                    if dt != np.inf:        # Check collision is possible.
                        self._pq[Event((j,self._N_balls,self._balls[j]._count,-1,dt + self._global_time))]\
                            = (dt + self._global_time)
                            
    def init_next_event(self):
        """ init_next_event | Initialisation of new event.
            Accounts for multiple simultaneous collisions. """
        self._events = []   # List of next events.
        min_event = self._pq.popitem()[0]
        self._min = min_event.dt()       # Find next event.
        self._events.append(min_event)      # Add back to list.
        
        while len(self._pq) != 0:       # Check for multiple collisions.
            if self._pq.peekitem()[0].dt().dt == self._min_dt:
                self._events.append(self._pq.popitem()[0])
            else: break
            
    def next_event(self):
        """ next_event | Selection of next collision event. """
        self._events = []
        min_event = self._pq.popitem()[0]
        
        while len(self._pq) != 0:       # Checks validity of the event.
            min_A, min_B = min_event.ball_A(), min_event.ball_B()
            if min_B == self._N_balls:  # Collision with container.
                if min_event.count_A() != self._balls[min_A]._count:
                    min_event = self._pq.popitem()[0]   # Find next event.
                else: break
            else:                       # Collision with other balls.
                if (min_event.count_A() != self._balls[min_A]._count 
                and min_event.count_B() != self._balls[min_B]._count):
                    min_event = self._pq.popitem()[0]   # Find next event.
                else: break
    
    def move_balls(self):
        """ move_balls | Move balls to next collision time. """
        for ball in self._balls:
            ball.move(self._min_dt - self._global_time)
            
    def collide_balls(self, pressure, test_pressure, brownian): # Remove test_press.
        """ collide_balls | Carry out ball collisions. Alter velocities.
            PARAMETERS
                pressure(boolean): Determine if pressure data is recorded.
                test_pressure(boolean): Determine if pressure data is recorded.
                brownian(boolean): Determine if Brownian motion data is recorded.
        """
        for event in self._events:
            ball_A, ball_B = event.ball_A(), event.ball_B()
            
            if ball_B == self._N_balls:         # Collision with container.
                self._balls[ball_A].collide(self._container)
                self._N_container_collisions += 1
                
                if pressure or test_pressure:
                    self._dp_container.append([np.linalg.norm(self._balls[ball_A]._dp), self._global_time])
            
            else:                               # Collision with other balls.
                self._balls[ball_A].collide(self._balls[ball_B])
            
            if brownian:
                if ball_A == 0: record = True
                
            else: record = False

        self._N_collisions += 1
        if brownian:
            if record: self._record_brownian()
    
    #%% Simulation Recording
    ### SIMULATION RECORDING METHODS
    # Gives simulation methods for recording data.

    def record_pressure(self): # Revisit, and alter.
        """ record_pressure | Record the pressure for every _(*)_ container collisions.
            RAISES
                IndexError: Insufficient number of collisions.
            RECORDS
            pd.DataFrame [pressure, t]
                CONTAINS
                -> pressure(float): Systemic pressure.
                -> t(float): Time.
        """
        self._pressure, N_collisions = [], 50   # Define (*)'s value.
        
        if not isinstance(self._dp_container, np.ndarray):
            self._dp_container = np.array(self._dp_container)
        
        try: max_t = self._dp_container[-1,1]
        except IndexError:      # No container collisions.
            print("Insufficient number of collisions for pressure calculation.")
            self._pressure = np.nan
            return
        
        min_t = max_t * 0.2     # Select final 20% of data at steady-state.
        start = 0
        while self._dp_container[start, 1] <= min_t:    # Find new start index.
            start += 1
            if start == len(self._dp_container - 1):
                print("Insufficient number of collisions.")
                self._pressure = np.nan
                return

        start += (len(self._dp_container) - start) % N_collisions
        dp_new = self._dp_container[start:,:]
        N_pressure = int(len(dp_new)/N_collisions)
        
        for i in range(N_pressure):
            index = i * N_collisions
            # <Pressure> = Σdp / T (period)
            pressure = np.sum(dp_new[index:index+N_collisions-1,0] /\
                (dp_new[index+N_collisions-1,1]-dp_new[index,1])*2*np.pi*self._r_container)    
            time = (dp_new[index+N_collisions-1,1] + dp_new[index,1]) / 2
            self._pressure.append([pressure,time])
        
        self._pressure = pd.DataFrame(self._pressure, columns = ['pressure', 't'])
        
    def record_pressure_moyen(self):
        """ record_pressure_moyen | Records the average systemic pressure.
            N.B: Must enable Simulation.run(pressure = True).
            RAISES
                IndexError: Insufficient number of collisions.
            RECORDS
                (float): average steady-state systemic pressure.
        """
        if not isinstance(self._dp_container, np.ndarray):
            self._dp_container = np.array(self._dp_container)
        
        try: max_t = self._dp_container[-1,1]
        except IndexError:      # No container collisions.
            print("Insufficient number of collisions for pressure calculation.")
            self._pressure = np.nan
            return
        
        min_t = max_t * 0.2     # Select final 20% of data at steady-state.
        start = 0
        while self._dp_container[start, 1] <= min_t:    # Find new start index.
            start += 1
            if start == len(self._dp_container - 1):
                print("Insufficient number of collisions.")
                self._pressure_moyen = np.nan
                return
        min_t = self._dp_container[start,1]
        # <Pressure> = Σdp / T (period)        
        self._pressure_moyen = np.sum(self._dp_container[start:,0]) /\
            ((max_t-min_t) * 2 * np.pi * self._r_container)

        
    def record_temperature(self):
        """ record_temperature | Records the temperature of the system, for all collisions.
            RECORDS
            -> pd.DataFrame [T, t, collision]
                CONTAINS
                -> T(float): Systemic temperature.
                -> t(float): Time.
                -> collision(int): Collision number.
        """
        KE = np.zeros(self._N_balls)
        for i, ball in enumerate(self._balls):
            KE[i] = 0.5 * ball._mass * np.linalg.norm(ball._vel)**2
        temperature = (np.sum(KE))/(self._N_balls*spc.Boltzmann)
        self._temperature.append([temperature, self._global_time, self._N_collisions])
        
        if len(self._KE) == self._collisions + 1:
            self._KE = pd.DataFrame(self._KE, columns=['T','t','collision'])
            
    def record_temperature_moyen(self):
        """ record_temperature_moyen | Records the average temperature of the system.
            RECORDS
                (float): Average systemic temperature.
        """
        self._temperature_moyen = np.mean(self._temperature['T'])

    def record_speed(self):
        """ record_speed | Record speed for all balls, for all collisions.
            RECORDS
                (list(float)): Speed of balls, for all collisions.
        """
        for ball in self._balls:
            self._speed.append(np.sqrt(np.dot(ball._vel,ball._vel)))

    def record_KE(self):
        """ record_KE | Record the total_KE of the system, for all collisions.
            RECORDS
            -> pd.DataFrame [KE, t, collision]
                CONTAINS
                -> KE(float): Systemic kinetic energy.
                -> t(float): Time.
                -> collision(int): Collision number.
        """    
        total_KE = np.sum([0.5 * self._m_balls * np.dot(ball._vel,ball._vel) for ball in self._balls])
        self._KE.append([total_KE, self._global_time, self._N_collisions])
        
        if len(self._KE) == self._collisions + 1:
            self._KE = pd.DataFrame(self._KE, columns=['KE','t','collision'])
        
    def record_distance_absolute(self):
        """ record_distance_absolute | Record the distance of all balls to O, ∀ collisions.
            RECORDS
                (list(float)): Distance of balls from O, for all collisions.
        """
        for ball in self._balls:
            self._distance_absolute.append(np.sqrt(np.dot(ball._pos,ball._pos)))

    def record_distance_relative(self):
        """ record_distance_relative | Record the distance between all balls, for all collisions.
            RECORDS
                (list(float)): Relative distance between ball pairs, for all collisions.
        """
        for _, pair in np.ndenumerate(self._pairs):
            ball_A, ball_B = pair[0], pair[1]
            distance_relative = np.sqrt(np.dot(ball_A._pos,ball_A._pos))\
                - np.sqrt(np.dot(ball_B._pos,ball_B._pos))
            self._distance_relative.append(distance_relative)
            
    def record_brownian(self, df=False):
        """ record_brownian | Records data prerequisite for Brownian motion investigation.
            RECORDS
            -> pd.DataFrame [x, y, t, collision]
                CONTAINS
                -> x(float): x co-ordinate of ball.
                -> y(float): y co-ordinate of ball.
                -> t(float): Collision time.
                -> collision(float): Collision number.
        """
        if df == False:
            self._brownian.append(np.array([self._ball[0]._pos[0],\
                self._ball[0]._pos[1], self._global_time, self._N_collisions]))
        else:       # Once data completes, add to pd.DataFrame
            self._brownian = pd.DataFrame(self._brownian,\
                columns = ['x','y','t','collision'])
            
    def record_dataset(self):
        """ record_dataset | Write simulation data into pd.DataFrame.
            RECORDS
            ->pd.DataFrame [ball, mass, x, y, v_x, v_y, collision, t, container]
                CONTAINS
                -> ball(int): Ball number.
                -> mass(float): Ball mass.
                -> x(float): Ball x co-ordinate.
                -> y(float): Ball y co-ordinate.
                -> v_x(float): Ball x velocity.
                -> v_y(float): Ball y velocity.
                -> collision(int): Collision number.
                -> t(float): Collision time.
                -> container(boolean): Detect if ball collides with container.
        """
        for i, ball in enumerate (self._balls):
            j = self._N_collisions * self._N_balls
            self._dataset[i+j, 0] = i
            self._dataset[i+j, 1] = ball._mass
            self._dataset[i+j, 2] = ball._pos[0]
            self._dataset[i+j, 3] = ball._pos[1]
            self._dataset[i+j, 4] = ball._vel[0]
            self._dataset[i+j, 5] = ball._vel[1]
            self._dataset[i+j, 6] = self._N_collisions
            self._dataset[i+j, 7] = self._global_time

            if (np.sqrt(np.dot(ball._pos,ball._pos)) - self._r_container + ball._radius <= 10e-10):       # Takes into account floating point.
                self._dataset[i+j, 8] = True
            else: self._dataset[i+j, 8] = False
        
        if self._N_collisions == self._collisions:
            self._dataset = pd.DataFrame(
                self._dataset,
                columns = ['ball','mass','x','y','v_x','v_y','collision',\
                    't','container']
            )
                
    def record_data_states(self,
                           distance_absolute = False, 
                           distance_relative = False, 
                           speed = False, KE = False,
                           temperature = False,
                           dataset = False):
        """ record_data_states | Records simulation datasets.
            PARAMETERS
                -> distance_absolute(boolean, optional): If True, writes dataset for all distances of ball to O, for all collisions.
                -> distance_relative(boolean, optional): If True, writes dataset for all relative distances between balls, for all collisions.
                -> speed(boolean, optional): If True, does the same as above for speeds of all balls, for all collisions.
                -> KE(boolean, optional): If True, same as above, for systemic KE.
                -> temperature(boolean, optional): If True, same as above, for the average temperature of the system (not for all collisions).
                -> dataset(boolean, optional): If True, same as above, for all simulation information.
        """      
        if distance_absolute:   self.record_distance_absolute()
        if distance_relative:   self.record_distance_relative()
        if speed:               self.record_speed()
        if KE:                  self.record_KE()
        if temperature:         self.record_temperature()
        if dataset:             self.record_dataset()

    def record_data_pressures(self, pressure=False, test_pressure=False):
        """ record_data_pressures | Records systemic pressure datasets.
            PARAMETERS
                -> pressure(boolean, optional): If True, writes dataset for the average systemic pressure.
                -> test_pressure(boolean, optional): If True, writes dataset of pressure for every _(*)_ container collisions.
        """ 
        if pressure:        self.record_pressure_moyen()
        if test_pressure:   self.record_pressure()

    #%% Simulation Property
    ### SIMULATION PROPERTY METHODS
    # Gives all the simulation methods about properties of the system.
    """
    Various properties and attributes of the system are returned.
    
    METHODS:
    N_balls | Returns how many container balls.
        RETURNS
            (int): Number of balls in container.
            
    balls | Lists all active balls in the simulation.
        RETURNS
            (list(bl.Ball)): a list of all the active balls in the simulation.
    
    container | Provides the container of the simulation.
        RETURNS
            (bl.Container): the simulation's container object.
    ________________________________________________________________________
    pressure | Provides the pressure values for every __ container collisions.
        N.B: Must enable Simulation.run(pressure = True).
        RETURNS
            pressure(float): system pressure.

    pressure_moyen | Provides the average pressure of the system.
        N.B: Must enable Simulation.run(pressure = True).
        RETURNS
            pressure_moyen(float): the system's average steady-state pressure.

    temperature | Provides the system temperature at all collision times.
        N.B: Must enable Simulation.run(temperature = True).
        RETURNS
        
    temperature_moyen | Provides the average temperature of the system.
        N.B: Must enable Simulation.run(temperature = True).
        RETURNS
            temperature_moyen(float): the system's average temperature.

    KE | Provides the total system KE, for all collisions.
        N.B: Must enable Simulation.run(KE = True).
        RETURNS
    
    distance_absolute | Provides distances from the origin for all balls, in all collisions.
        N.B: Must enable Simulation.run(distance_absolute = True).
        RETURNS
            distance_absolute(list(float)): ball distances from the origin, for all collisions.

    distance_relative | Provides relative distance between all balls, in all collisions.
        N.B: Must enable Simulation.run(distance_relative = True).
        RETURNS
            relative(list(float)): ball distances between all possible pairs of balls, in all collisions.
            
    speed | Provides the speed for all balls, in all collisions.
        N.B: Must enable Simulation.run(speed = True).
        RETURNS
            speed(list(float)): ball speeds, in all collisions.
            
    brownian | Provides the Brownian motion investigation dataset.
        N.B: Must enable Simulation.run(brownian = True).
        RETURNS
    
    dataset | Provides a complete dataset, for this simulation.
        N.B: Must enable Simulation.run(dataset = True).
        RETURNS
            See record_dataset() for more information.
        
    pair_combn | Provides a complete list of all possible combinations of ball pairs.
        PARAMETERS
            container (boolean, optional): include Container in pairs.
        RETURNS
            (list(tuple(int))): list containing all tuples of pairs.
    """
    def N_balls(self):              return self._N_balls
    def balls(self):                return self._balls
    def container(self):            return self._container
    
    def pressure(self):             return self._pressure
    def pressure_moyen(self):       return self._pressure_moyen
    def temperature(self):          return self._temperature
    def temperature_moyen(self):    return self._temperature_moyen
    def KE(self):                   return self._KE
    def distance_absolute(self):    return self._distance_absolute
    def distance_relative(self):    return self._distance_relative
    def speed(self):                return self._speed
    def brownian(self):             return self._brownian
    def dataset(self):              return self._dataset
    def pair_combn(self, container=False):
        if not container:       list_number = list(range(self._N_balls))
        else:                   list_number = list(range(self._N_balls + 1))
        return                  list(it.combinations(list_number, 2))
   
    #%% Simulation Randomisation
    ### SIMULATION RANDOMISATION METHODS
    # Gives all the randomisation methods for the simulation.
    def gen_random_uniform(maximum_range):
        """ gen_random_uniform | Provides a uniform distribution centered at (0,0), generating random floats.
            PARAMETERS
                maximum_range (float): sets the maximum range for the distribution.
            RETURNS
                (float): a random float between [-max_range, max_range].
        """
        return np.random.uniform(-maximum_range,maximum_range)   

    def gen_random_velocities(self, N_balls, maximum_speed_range):
        """ gen_random_velocities | Generates random velocities from a uniform distribution
        for a given number of balls in the range of [-maximum_speed_range, maximum_speed_range].
            PARAMETERS
                N_balls (int): the number of balls.
                random_speed_range(float): the ± range in the x, and y-components.
            RETURNS
        """
        list = []
        for _ in range(N_balls):
            list.append(np.array(self.gen_random_positions(maximum_speed_range),\
                self.gen_random_uniform(maximum_speed_range)))
        return list

    def gen_random_positions(self, start=0):
        """ gen_random_positions | Generates the non-overlapping random ball positions.
            PARAMETERS
                start (boolean, optional):
            RAISES
                Exception:
                    Balls cannot fit in this container. Reduce N_balls / increase r_container.
        """
        for i in range(start, self._N_balls):
            position, error_count = np.zeros(2), 0
            while True:
                # Extreme case error handling to prevent computational overload.
                if error_count > 1e5:
                    raise Exception \
                        ("The area of this container is too small for ball size.")
                x = self.gen_random_uniform(self._r_container - self._balls[i]._radius)
                y = self.gen_random_uniform(self._r_container - self._balls[i].radius)
                # Check if the randomly-assigned position is valid.
                while(np.sqrt(x**2 + y**2) >= self._r_container - self._balls[i]._radius):
                    x = self.gen_random_uniform(self._r_container - self._balls[i]._radius)
                    y = self.gen_random_uniform(self._r_container - self._balls[i].radius)
                position = np.array([x,y])
                append = False
                
                for j in range(0, i):
                    distance = np.sqrt((self._balls[j]._pos[0] - position[0])**2 \
                        + (self._balls[j]._pos[1] - position[1]**2))
                    if distance <= self._balls[i].radius + self._balls[j]._radius:
                        append = False
                        error_count += 1
                        break
                    else: append = True
                if append or i == 0: break
            self.ball[i].set_pos(position)

    #%% Simulation Brownian
    ### SIMULATION BROWNIAN MOTION INVESTIGATION
    # Gives all the methods for the Brownian Motion investigation.
    
    def brownian_init(self, radius = 5, mass = 10):
        """ brownian_init | Initialisation of Brownian motion investigation, ball "0".
            -> Position of ball 0 is initialised to (0,0).
            -> Remaining balls are randomly distributed in the container.
            PARAMETERS
                radius (float, optional): radius of ball being investigated.
                mass (float, optional): mass of ball being investigated.
        """
        self._balls[0].set_pos(np.array([0.0,0.0]))
        self._balls[0].set_radius(radius)
        self._balls[0].set_mass(mass)
        self.gen_random_positions(start = 1)

    def set_ball_velocities(self, velocity_list):
        """ set_ball_velocities | Sets the velocity of all the balls from a list of velocities.
        PARAMETERS
            velocity_list (list of (np.ndarray of (floats))): lists all the ball velocities in their x- and y- directions.
        """
        for i, velocity in np.ndenumerate(velocity_list):
            self._balls[i].set_vel(velocity)

    def brownian_velocities(self, maximum_speed):
        """ brownian_velocities | Generates and assigns random velocities to all balls from a uniform random 
        distribution of x, and y-components of velocity.
            PARAMETERS
                maximum_speed (float): the range to generate velocities from.
        """                    
        list = self.gen_random_velocities(self._N_balls, self._random_speed)
        self.set_ball_velocities(list)
        
    def brownian_tracer(self):
        """ brownian_tracer | Draws out path in animation followed by ball "0".
            RETURNS
                (plt.Line2D): path travelled between collisions.
        """
        trace = plt.Line2D(
            xdata=[self._balls[0]._pos[0], self._brownian[-1][0]],
            ydata=[self._balls[0]._pos[1], self._brownian[-1][1]],
            color="black", alpha=0.6,lw=0.7)
        return trace

    #%% Simulation Display
    ### SIMULATION DISPLAY METHODS
    # Gives all the links with Pylab and Matplotlib.

    def patch_init(self):
        """ patch_init | Initialisation of the ball and container patches as part of the animation. 
        The balls and container are drawn using plt.pyplot.Circle objects.
        """
        ball_patches = [] # Creates a list of containing all ball patches.
        position_container = self._container._pos
        r_container = self._r_container
        container_outline = plt.Circle(position_container,r_container, ec="b", fill=False, ls="solid")
        
        for i, ball in np.ndenumerate(self._balls): ## Change to enumerate() fail.
            pos, radius_ball = ball._pos, ball._radius
            
            if i != 0:
                ball_patches.append(\
                    plt.Circle(pos, radius_ball, ec="black", fc=tuple((np.random.rand(),np.random.rand(),np.random.rand()))))
            else: # Set the first ball to yellow - change later.
                ball_patches.append(\
                    plt.Circle(pos, radius_ball, ec="black", fc=tuple(np.random.rand(),np.random.rand(),np.random.rand())))
        self._ball_patches = ball_patches
        self._container_outline = container_outline
    
    def draw(self):
        """ draw | Draw the current static simulation state.
        """
        self.patch_init()
        
        plt.figure(num="Current Simulation State")
        ax = plt.axes(
            xlim=(-self._r_container, self._r_container),
            ylim=(-self._r_container,self._r_container),
            aspect="equal",)
        ax.add_patch(self._container_outline)
        for patch in self._ball_patches:
            ax.add_patch(patch)
        plt.show()
        
    def update_patch(self):
        """ update_patch | Updates the animation with new positions of ball patches.
        """
        for i in range(0,self._N_balls):
            self._ball_patches[i].center = self._balls[i].pos()
    
    #%% Simulation Miscellaneous
    ### SIMULATION MISCELLANEOUS METHODS
    def linear_inter(x, x_1, x_2, y_1, y_2):
        """ linear_inter | Returns y-value for two given linear equations.
        """
        return y_1 + (y_2 - y_1) / (x_2 - x_1) * (x - x_1)
      
    def brownian_equal_samples(self, data, N_samples):
        """ brownian_equal_samples | Samples positions of the ball given regular time intervals.
        PARAMETERS
            data: pd.DataFrame [x,y,t]
            -> x(float): Ball's x co-ordinate.
            -> y(float): Ball's y co-ordinate.
            -> t(float): Collision time.
            N_samples(int): Number of samples from data.
        RETURNS
            pd.DataFrame [x,y,t]
            -> x(float): Ball's x co-ordinate.
            -> y(float): Ball's y co-ordinate.
            -> t(float): Equally-sampled collision times.
        """
        x_pos, x_samp = np.array(data["x"]), np.zeros(N_samples)
        y_pos, y_samp = np.array(data["y"]), np.zeros(N_samples)
        time, t_samp = np.array(data["t"]), np.zeros(N_samples)

        dt = (time[-1] - time[0])/(N_samples + 1)
        index, i = 0, 0
        
        # Sample time positions for interval dt. Use np.linspace
        while i < N_samples:
            if (time + i * dt) == time[index]:
                x_samp[i], y_samp[i] = x_pos[index], y_pos[index]
                t_samp[i] = time[index]; i += 1
            else:       # Time exceeds next indexed value.
                x1, x2 = x_pos[index], x_pos[index + 1]
                y1, y2 = y_pos[index], y_pos[index + 1]
                t1, t2 = time[index], time[index + 1]
                
                x, y = self.linear_inter(t_samp,t1,t2,x1,x2), self.linear_inter(t_samp,t1,t2,y1,y2)
                x_samp[i], y_samp[i], t_samp[i] = x, y, t_samp
                i += 1
                
        return pd.DataFrame([x_samp, y_samp, t_samp]).transpose().rename(columns=\
            {0: "x", 1: "y", 2: "t"})
        
    def brownian_paths_tracer(data):
        """ brownian_paths_tracer | Trace out the path(s) of a given ball with positional data.
        PARAMETERS
            data: pd.DataFrame [x,y] - the time-ordered positional data at collision.
            -> x(float): Ball's x co-ordinate.
            -> y(float): Ball's y co-ordinate.
        RETURNS
            list(plt.Line2D): List of plottable Line2D objects.
        """
        x, y, path_list = ["x"], ["y"], []
        for i in range(1, len(x)):
            path = plt.Line2D(xdata=[x[i],x[i-1]], ydata=[y[i],y[i-1]],
                              color='0.01', alpha=0.02)
            path_list.append(path)
        return path_list

    titleFont =     {'fontname': 'Kinnari', 'size': 13}
    axesFont =      {'fontname': 'Kinnari', 'size': 9}
    ticksFont =     {'fontname': 'SF Mono', 'size': 7}
    errorStyle =    {'mew': 1, 'ms': 3, 'capsize': 3, 'color': 'blue', 'ls': ''}
    pointStyle =    {'mew': 1, 'ms': 3, 'color': 'blue'}
    lineStyle =     {'linewidth': 0.5}
    lineStyleBold = {'linewidth': 1}
    histStyle =     {'facecolor': 'green', 'alpha': 0.5, 'edgecolor': 'black'}
    
    #%% Simulation Run
    ### SIMULATION RUN METHOD
    # The method to run simulations. 
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
        distance_absolute = False,  # Use ball distance from the origin.
        distance_relative = False,  # Use relative distances between balls.
        progress_bar = True,        # Enable progress bar animation in terminal.
        animate = False             # Enables animation.
    ):
        """ run | Runs the 2D simulation of colliding particles within the container.
        PARAMETERS
        -> collisions (int, optional): number of collisions in the simulation.
        -> time (float, optional): time period (s) between animation frames.
            
        -> pressure (boolean, optional): records pressure for every __ collisions with the wall of the container.
        -> temperature (boolean, optional): records temperature.
        -> KE (boolean, optional): records the system KE for every collision.
        -> speed (boolean, optional): records speed of all balls in all collisions.
        -> brownian (boolean, optional): records data for Brownian motion.
        -> dataset (boolean, optional): records dataset of simulation information.
        -> distance_absolute (boolean, optional): records the distances of all balls from the origin, in all collisions.
        -> distance_relative (boolean, optional): records the relative distances between all the balls, in all collisions.
        -> progress_bar (boolean, optional): displays a progress bar.

        RETURNS
            gloss(dict): Glossary of needed datasets.
        """
        self._global_time = 0
        self._collisions = collisions
        self._pq = hd.heapdict()
        
        print("Starting {self._N_balls} balls, r_balls = {self._r_balls}, \
            speed range = {self._speed_range}, r_container = {self._r_container}, \
            {self._collisions} collisions.")
        
        if animate:                 # Initialise animation.
            self.patch_init()
            plt.figure(num="Simulation Animation")
            plt.rcParams.update(plt.rcParamsDefault)
            ax = plt.axes(xlim=(-self._r_container,self._r_container),
                          ylim=(-self._r_container,self._r_container),
                          aspect="equal") # Equal square axes of r_container.
            
            ax.add_patch(self._container_outline)
            for ball_patch in self._ball_patches: ax.add_patch(ball_patch)
            plt.pause(time)

        if dataset: self._dataset = np.zeros((self._N_balls*(self._collisions*1),9))
        if brownian: self.record_brownian()
        
        self.record_data_states(distance_absolute=distance_absolute,
                                distance_relative=distance_relative,
                                speed=speed, KE=KE,
                                temperature=temperature,
                                dataset=dataset)

        self.init_collision_time()      # Run the first collision.
        self.init_next_event()
        self.move_balls()
        self._global_time = self._min_dt
        
        if animate:
            self.update_patch()
            if brownian: ax.add_line(self.brownian_tracer())
            plt.pause(time)
        
        self.collide_balls(pressure,test_pressure,brownian)
        self.record_data_states(distance_absolute=distance_absolute,
                                distance_relative=distance_relative,
                                speed=speed, KE=KE,
                                temperature=temperature,
                                dataset=dataset)
        
        if progress_bar: self._time_epoch = tm.time()

        for i in range(2,collisions+1):
            if progress_bar: progress_bar(self._time_epoch,i,collisions)
            self.collision_time()
            self.next_event()
            self.move_balls()
            
            self._global_time = self._min_dt
            
            if animate: 
                self.update_patch()
                if brownian: ax.add_line(self.brownian_tracer())
                plt.pause(time)
        self.collide_balls(pressure,test_pressure,brownian)

        self.record_data_states(distance_absolute=distance_absolute,
                                distance_relative=distance_relative,
                                speed=speed, KE=KE,
                                temperature=temperature,
                                dataset=dataset)
        
        if animate: plt.show()
        if temperature: self.record_temperature_moyen()
        
        self.record_data_pressures(pressure=pressure, test_pressure=test_pressure)
    
        if brownian: self.record_brownian(df = True)        

        gloss = self.glossary(speed=speed, KE=KE,
                              distance_absolute=distance_absolute,
                              distance_relative=distance_relative,
                              temperature=temperature,
                              test_temperature=test_temperature,
                              pressure=pressure,
                              test_pressure=test_pressure,
                              dataset=dataset,
                              brownian=brownian
                              )
        
        print("Ending {self._N_balls} balls, r_balls = {self._r_balls}, \
            speed range = {self._speed_range}, r_container = {self._r_container}, \
            {self._collisions} collisions.")
        return gloss
        
class Event(tuple):
    """
        EVENT CLASS | 
        A tuple of 5 elements (ball_A, ball_B, count_A, count_B, dt).
        PARAMETERS
            ball_A (int): The first ball in impending collision.
            ball_B (int): The second ball in impending collision.
            count_A (int): The number of collisions the first ball 
                encountered prior to this impending collision calculation.
            count_B (int): The number of collisions the second ball 
                encountered prior to this impending collision calculation.
            dt (float): The global time this collision will happen on. 
        RETURNS
            ball_A (int): the index of first ball in impending collision.
            ball_B (int): the index of second ball in impending collision.
            count_A (int): the number of collisions the first ball encountered prior to this impending collision calculation.
            count_B (int): the number of collisions the second ball encountered prior to this impending collision calculation.
            dt (float): the global time(step) this collision will happen on.
            pair (list(int)): a list of the two balls or a ball and container involved in the collision.
    """
    def ball_A(self):   return self[0]
    def ball_B(self):   return self[1]
    def count_A(self):  return self[2]
    def count_B(self):  return self[3]
    def dt(self):       return self[4]
    def pair(self):     return [self[0], self[1]]
    
    