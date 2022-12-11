""" SIMULATION MODULE (import simulation as sm) | Created by 𝑀𝒶𝓇𝓉𝒾𝓃 𝒜. 𝐻𝑒, 2022.11.20
    For the simulation of elastic collisions of balls with others, and container.
    -> Simulation Class.
    -> Event Class.
    -> Plots Class.
""" # martinhe.com/thermo-project
import ball as bl
import numpy as np
import scipy as sp
import scipy.optimize as opt
import scipy.constants as spc
import matplotlib.pyplot as plt
import matplotlib.font_manager as fnt
import seaborn as sns
import pandas as pd
import itertools as it
import heapdict as hd
import time as tm
import concurrent.futures as ft
import sys, os, random, warnings

titleFont =     {'fontname': 'C059', 'size': 13}
axesFont =      {'fontname': 'C059', 'size': 9}
ticksFont =     {'fontname': 'SF Mono', 'size': 7}
errorStyle =    {'mew': 1, 'ms': 3, 'capsize': 3, 'color': 'blue', 'ls': ''}
pointStyle =    {'mew': 1, 'ms': 3, 'color': 'blue'}
lineStyle =     {'linewidth': 0.5}
lineStyleBold = {'linewidth': 1}
histStyle =     {'facecolor': 'green', 'alpha': 0.5, 'edgecolor': 'black'}

def gen_random_uniform(maximum_range):      # Please keep global.
    """ gen_random_uniform | Provides a uniform distribution centered at (0,0),
        generating random floats.
            <PARAMETERS>
                maximum_range (float): sets the maximum range for the distribution.
            RETURNS
                (float): a random float between [-max_range, max_range].
    """
    return np.random.uniform(-maximum_range,maximum_range)

class Simulation:
    """ SIMULATION CLASS |
        Simulate movement of hard spherical gas particles in circular container.
            < PARAMETERS >
            -> N_balls(int, optional): The number of balls to simulate.
            -> r_balls(float, optional): The radius of balls in the simulation.
            -> m_balls(float, optional): The mass of balls in the simulation.
            -> r_container(float, optional): The radius of the circular container.
            -> random_position(boolean, optional): Randomly-position balls.
            -> random_speed_range(float, optional): Range of speeds ball can take.
            < VARIABLES >
                -> balls(list[bl.Ball]): The list containing all Ball objects.
                -> N_balls(int): The number of simulation balls.
                -> N_collisions(int): The collision count incrementer.
                -> m_balls(float): The mass of simulation balls.
                -> r_balls(float): The radius of simulation balls.
                -> r_container(float): The container radius.
                -> pairs(list([bl.Ball,bl.Ball])): List of Ball pair tuples.
                -> temperature(list(float)): List of system temperature over time.
                -> KE(list(float)): List of system kinetic energy over time.
                -> speed(list(float)): List of speed of Ball objects in Simulation.
                -> distance_absolute(list(float)): List of all Ball object distance 
                                                   from origin, i.e. the centre.
                -> distance_relative(list(float)): List of distances between all
                                                   pairs of Ball objects.
                -> random_position(boolean): Enable random positioning of Ball.
                -> random_speed_range(int): Set range of random speeds of Ball.
                -> brownian(boolean): Enables the Brownian pollen investigation.
                -> container(bl.Container): The Container object for Simulation.
                -> dp_container(list(float)): List of container momentum changes.
                -> N_container_collisions(int): Number of container collisions.
                -> utc(float): Universal Co-ordinated Time for Simulation.
                -> min_dt(Event): Event with minimmum next time to collision.
                -> pq(hd.heapdict): Dictionary containing priority queue.
                -> events(list[Event]): List of current Event objects.
    """
    def __init__(self,
    # 𝔅𝔞𝔩𝔩 𝔓𝔯𝔬𝔭𝔢𝔯𝔱𝔦𝔢𝔰:
        N_balls = 1,                        # The number of balls.
        r_balls = 1,                        # The radius of the balls.
        m_balls = 1,                        # The mass of the balls.
        r_container = 10,                   # The radius of the container.
    # 𝔅𝔞𝔩𝔩 ℜ𝔞𝔫𝔡𝔬𝔪𝔦𝔰𝔞𝔱𝔦𝔬𝔫:
        random_position = True,             # Randomly-positioned balls.
        random_speed_range = 5              # Range for random speed generation.
    ):
    # 𝔅𝔞𝔩𝔩 𝔓𝔯𝔬𝔭𝔢𝔯𝔱𝔦𝔢𝔰:
        self._balls = []                    # List of container balls.
        self._N_balls = N_balls             # The number of container balls.
        self._N_collisions = 0              # The collision count incrementer.
        self._m_balls = m_balls             # The mass of the container balls.
        self._r_balls = r_balls             # The radius of the container balls.
        self._r_container = r_container     # The radius of the container.
        self._pairs = self.property(pairs=True)     # List of ball pairs.
    # ℭ𝔬𝔩𝔩𝔦𝔰𝔦𝔬𝔫 𝔓𝔯𝔬𝔭𝔢𝔯𝔱𝔦𝔢𝔰:
        self._temperature = []              # System T, ∀ collisions.
        self._KE = []                       # System K.E, ∀ collisions.
        self._speed = []                    # Speed of all balls ∀ collisions.
        self._distance_absolute = []        # Ball distance from origin.
        self._distance_relative = []        # Relative distance between balls.
    # 𝔅𝔞𝔩𝔩 ℜ𝔞𝔫𝔡𝔬𝔪𝔦𝔰𝔞𝔱𝔦𝔬𝔫:
        self._random_position = random_position     # Set random positioning.
        self._random_speed = random_speed_range     # Set random speed range.
        self._brownian = []                 # Brownian investigation dataset.
    # ℭ𝔬𝔫𝔱𝔞𝔦𝔫𝔢𝔯 𝔓𝔯𝔬𝔭𝔢𝔯𝔱𝔦𝔢𝔰:
        self._container = bl.Container(radius=r_container) # Chosen container.
        self._dp_container = []             # Changes in container momentum.
        self._N_container_collisions = 0    # Number of container collisions.
    # 𝔗𝔦𝔪𝔢 𝔓𝔯𝔬𝔭𝔢𝔯𝔱𝔦𝔢𝔰:
        self._utc = 0                       # Universal Co-ordinated Time.
        self._min_dt = 0                    # Minimum time to next collision.
    # 𝔈𝔳𝔢𝔫𝔱 𝔓𝔯𝔬𝔭𝔢𝔯𝔱𝔦𝔢𝔰:
        self._queue = hd.heapdict()         # Dictionary of queued objects.
        self._events = []                   # List of current events.

        """ __init__ | The intialisation method for Simulation.
        Creates a new Simulation object. For properties, see Simulation.
        """
        for _ in range(0, N_balls):         self._balls.append(bl.Ball(radius=\
                                                r_balls, mass=m_balls))
        if random_position:                 self.gen_random_positions()
        self.gen_random_velocities(max_speed=random_speed_range, N_balls=N_balls)
        warnings.filterwarnings('ignore')

    def __repr__(self):
        """ __repr__ | Gives the representation format for Simulation.
            RETURNS
                __repr__(string): Contains Simulation properties in __repr__.
        """
        return ("SIMULATION PROPERTIES: N = %s balls, m_balls = %s,\
                r_balls = %s, r_container = %s")%(self._N_balls, self._m_balls,\
                    self._r_balls, self._r_container)

    def __str__(self):
        """ __str__ | Gives the string format for Simulation.
                RETURNS
                    __str__(string): Contains Simulation properties in __str__.
        """
        return ("SIMULATION PROPERTIES: N = %s balls, m_balls = %s,\
                r_balls = %s, r_container = %s")%(self._N_balls, self._m_balls,\
                    self._r_balls, self._r_container)

    def glossary(self, pressure=False, test_pressure=False, temperature=False,
                 test_temperature=False, speed=False, KE=False,
                 distance_absolute=False, distance_relative = False,
                 dataset=False, brownian=False):
        """ glossary | Append glossary of various datasets at end of simulation.
                < PARAMETERS >
                All parameters are self-explanatory booleans. For more info,
                see Documentation in __init__ method, or record_property function.
        """
        gloss = {}
        if speed: gloss["speed"] = self._speed
        if KE: gloss["KE"] = self._KE
        if distance_absolute: gloss["distance from centre"] = self._distance_absolute
        if distance_relative: gloss["relative distance"] = self._distance_relative
        if temperature: gloss["average temperature"] = self._temperature_moyen
        if test_temperature: gloss["temperature"] = self._temperature
        if pressure: gloss["average pressure"] = self._pressure_moyen
        if test_pressure: gloss["pressure"] = self._pressure
        if dataset: gloss["dataset"] = self._dataset
        if brownian: gloss["brownian"] = self.brownian
        return gloss
  
    def property(self, N_balls=False, ball=False, container=False, speed=False,
                 KE=False, distance_absolute=False, distance_relative=False,
                 temperature=False, temperature_moyen=False, pressure=False,
                 pressure_moyen=False, brownian=False, dataset=False, pairs=False,
                 pairs_container=False):
        """ property | Return various properties of the simulation.
            < PARAMETERS >
            -> N_balls: Returns how many container balls.
                RETURNS
                    (int): Number of balls in container.
            -> ball: Lists all active balls in the simulation.
                RETURNS
                    (list(bl.Ball)): a list of all the active balls in the 
                                     simulation.
            -> container: Provides the container of the simulation.
                RETURNS
                    (bl.Container): the simulation's container object.
            -> speed: Provides the speed for all balls, in all collisions.
                N.B: Must enable Simulation.run(speed = True).
                RETURNS
                    speed(list(float)): ball speeds, in all collisions.
            -> KE: Provides the total system KE, for all collisions.
                N.B: Must enable Simulation.run(KE = True).
                RETURNS
                    KE(list(float)): Total systemic KE, at all collisions.
            -> distance_absolute: Provides distances from the origin for all 
                                  balls, in all collisions.
                N.B: Must enable Simulation.run(distance_absolute = True).
                RETURNS
                    distance_absolute(list(float)): ball distances from the 
                    origin, for all collisions.
            -> distance_relative: Provides relative distance between all balls,
                                  in all collisions.
                N.B: Must enable Simulation.run(distance_relative = True).
                RETURNS
                    relative(list(float)): ball distances between all possible 
                                           pairs of balls, in all collisions.
            -> temperature: Provides the system temperature at all collision
                            times.
                N.B: Must enable Simulation.run(temperature = True).
                RETURNS
                    temperature(list(float)): Systemic temperature, at all 
                                              collision times.
            -> temperature_moyen: Provides the average temperature of the system.
                N.B: Must enable Simulation.run(temperature = True).
                RETURNS
                    temperature_moyen(float): the system's average temperature.
            -> pressure: Provides the pressure values for every __ container 
                         collisions.
                N.B: Must enable Simulation.run(pressure = True).
                RETURNS
                    pressure(float): system pressure.
            -> pressure_moyen: Provides the average pressure of the system.
                N.B: Must enable Simulation.run(pressure = True).
                RETURNS
                    pressure_moyen(float): the system's average steady-state 
                                           pressure.
            -> brownian: Provides the Brownian motion investigation dataset.
                N.B: Must enable Simulation.run(brownian = True).
                RETURNS
            -> dataset: Provides a complete dataset, for this simulation.
                    N.B: Must enable Simulation.run(dataset = True).
                    RETURNS
                        See record_dataset() for more information.
            -> pair_combn: Provides a complete list of all possible combinations
                           of ball pairs.
                    RETURNS
                        (list(tuple(int))): list containing all tuples of pairs.              
            -> pairs_container(boolean, optional): include Container in pairs.
        """

        if N_balls:              return self._N_balls              
        if ball:                 return self._balls
        if container:            return self._container
        if speed:                return self._speed
        if KE:                   return self._KE
        if distance_absolute:    return self._distance_absolute
        if distance_relative:    return self._distance_relative
        if temperature:          return self._temperature
        if temperature_moyen:    return self._temperature_moyen
        if pressure:             return self._pressure
        if pressure_moyen:       return self._pressure_moyen
        if brownian:             return self._brownian
        if dataset:              return self._dataset
        if pairs:
            if not pairs_container: list_number = list(range(self._N_balls))
            else:                   list_number = list(range(self._N_balls + 1))
            return                  list(it.combinations(list_number, 2))

    def set_vel_ball(self, velocity_list):
        """ set_ball_velocities | Sets the velocity of all the balls from a list
                                  of velocities.
                <PARAMETERS>
                    velocity_list (list of (np.ndarray of (floats))): lists all 
                    the ball velocities in their x- and y- directions.
        """
        for i, vel in enumerate(velocity_list):     self._balls[i].set_vel(vel)
        
    def next_collision(self, init=False, pressure=False, test_pressure=False, 
                       brownian=False):
        """ next_collision | Sets up and performs the next collision.
            The logic of the next_collision method is given by:
            1. Find the time to the next collision.
            2. Move the system to that point in time.
            3. Perform the collision.
                < PARAMETERS >
                -> init(boolean, optional): For the first collision.
                -> pressure(boolean, optional): Enables pressure calculations.
                -> test_pressure(boolean, optional): Enables test pressure.
                -> brownian(boolean, optional): Enables Brownian investigation.
        """
        if init:
            for pair in self._pairs:    # ∀ possible combinations of ball pairs.
                ball_A, ball_B = self._balls[pair[0]], self._balls[pair[1]]
                dt = ball_A.time_to_collision(ball_B)
                if dt != np.inf:        # Only consider possible collisions.
                    self._queue[Event((pair[0], pair[1], ball_A._count, \
                                    ball_B._count,dt))] = dt

            for i, ball in enumerate(self._balls):  # Ball-Container collisions.
                dt = ball.time_to_collision(self._container)
                if dt != np.inf:
                    self._queue[Event((i,self._N_balls,ball._count,-1,dt))] = dt
                    
            self._events = []                   # Clear list of next events.
            min_event = self._queue.popitem()[0]   # Select next event.
            self._min = min_event.dt()          # Find next event.
            self._events.append(min_event)      # Add back to list.

            while len(self._queue) != 0:           # Check for multiple collisions.
                if self._queue.peekitem()[0].dt() == self._min_dt:
                    self._events.append(self._queue.popitem()[0])
                else: break
        
        elif init==False:      
            collided_ball = set()
            for event in self._events:          # Add events for next collision(s).
                for collided in event.pair():
                    collided_ball.add(collided)

            for element in collided_ball:           # Add events to the heapdict.
                if element != self._N_balls:        # Collisions with the container.
                    dt = self._balls[element].time_to_collision(self._container)
                    if dt != np.inf:
                        self._queue[Event((element, self._N_balls, self._balls\
                                [element]._count, -1, dt+self._utc,))] =\
                                (dt + self._utc)

                    for j in range(self._N_balls):  # Collisions with other balls.
                        if j != element:
                            if j < element:         # Ensure correct index order.
                                ball_A, ball_B = self._balls[j],self._balls[element]
                                index_A, index_B = j, element
                            else:
                                ball_A, ball_B = self._balls[element],self._balls[j]
                                index_A, index_B = element, j

                            dt = ball_A.time_to_collision(ball_B)
                            if dt != np.inf:
                                self._queue[Event((index_A, index_B, self._balls\
                                    [index_A]._count, self._balls[index_B]._count,\
                                    dt + self._utc))] = (dt + self._utc)

                else:                               # Collision with container.
                    for j in range(self._N_balls):
                        dt = self._balls[j].time_to_collision(self._container)
                        if dt != np.inf:
                            self._queue[Event((j,self._N_balls,self._balls[j]._count,\
                                    -1, dt + self._utc))] = (dt + self._utc)

            self._events = []
            min_event = self._queue.popitem()[0]

            while len(self._queue) != 0:           # Check if event is possible.
                min_A, min_B = min_event.ball_A(), min_event.ball_B()
                if min_B == self._N_balls:      # Container collision.
                    if min_event.count_A() != self._balls[min_A]._count:
                        min_event = self._queue.popitem()[0]
                    else: break
                else:                           # Ball-Ball collision.
                    if (min_event.count_A() != self._balls[min_A]._count
                    and min_event.count_B() != self._balls[min_B]._count):
                        min_event = self._queue.popitem()[0]
                    else: break

            self._min_dt = min_event.dt(); self._events.append(min_event)

            while len(self._queue) != 0:           # Checks for multiple collisions.
                next_event = self._queue.peekitem()[0]
                if next_event.dt() == self._min_dt:
                    next_A, next_B = next_event.ball_A(), next_event.ball_B()
                    if next_B == self._N_balls: # Container collision.
                        if next_event.count_A() == self._balls[next_A]._count:
                            self._events.append(self._queue.popitem()[0])
                        else: break
                    else:                       # Ball-Ball collision.
                        if (next_event.count_A() == self._balls[next_A]._count
                            and next_event.count_B() == self._balls[next_B]._count):
                            self._events.append(self._queue.popitem()[0])
                        else: break
                else: break
                
            for ball in self._balls:     ball.move(self._min_dt - self._utc)
        
        self._utc = self._min_dt
        if self._animate:
            self.draw(init=False)
            if self._has_brownian: path = self.brownian_investigation(); self._axes.add_line(path)
            plt.pause(self._time)

        record = False
        for event in self._events:
            ball_1, ball_2 = event.ball_A(), event.ball_B()
            if ball_2 == self._N_balls:
                self._balls[ball_1].collide(self._container)
                self._N_container_collisions += 1
                if pressure or test_pressure:       # Append to dp_container.
                    self._dp_container.append([np.linalg.norm(self._balls\
                                              [ball_1]._dp), self._utc])
            else:  self._balls[ball_1].collide(self._balls[ball_2])
            if self._has_brownian:
                if ball_1 == 0: record = True

        self._N_collisions += 1
        if self._has_brownian:
            if record:  self.record_property(brownian=True)

    def gen_random_velocities(self, N_balls=50, max_speed=5):
        """ gen_random_velocities | Generate random velocities from a uniform 
            distribution for a given number of balls in the range of:
            [-maximum_speed_range, maximum_speed_range].
                <PARAMETERS>
                -> random_speed_range(float): the ± range in the x, and 
                                              y-components.
        """
        list = []
        for _ in range(N_balls):
            list.append(np.array([gen_random_uniform(max_speed),
                                  gen_random_uniform(max_speed)]))
        self.set_vel_ball(list)

    def gen_random_positions(self, start=0):
        """ gen_random_positions | Generate the non-overlapping random ball 
                                   positions.
                < PARAMETERS >
                -> start (float, optional): The starting position.
                RAISES
                    Exception: Balls cannot fit in this container. Reduce
                               N_balls or increase r_container.
        """
        for i in range(start, self._N_balls):
            position, error_count = np.zeros(2), 0
            while True:             # Error handling to prevent computer death.
                if error_count > 1e5:
                    raise Exception("The area of this container is too small\
                                     for ball size.")
                x = gen_random_uniform(self._r_container - \
                                       self._balls[i]._radius)
                y = gen_random_uniform(self._r_container - \
                                       self._balls[i]._radius)

                while(np.sqrt(x**2 + y**2) >= self._r_container - \
                      self._balls[i]._radius):  # Check if inside container.
                    x = gen_random_uniform(self._r_container - \
                                           self._balls[i]._radius)
                    y = gen_random_uniform(self._r_container - \
                                           self._balls[i]._radius)
                position, append = np.array([x, y]), False
                
                for j in range(0, i):           # Final check!
                    distance = np.sqrt((self._balls[j]._pos[0] - position[0])**2 \
                        + (self._balls[j]._pos[1] - position[1]) ** 2)
                    if distance <= self._balls[i]._radius + self._balls[j]\
                        ._radius: append = False; error_count += 1; break
                    else:               append = True
                if append or i == 0:    break
            self._balls[i].set_pos(position)

    def draw(self, init=False, draw_current_state=False, random_color=True):
        """ draw | Updates the animation with new positions of ball patches.
                < PARAMETERS >
                -> init(boolean, optional): Draw first Simulation state.
                -> draw_current_state(boolean, optional): Draw current state.
        """
        def init_patches(self):     self.draw(init=True)
        if init:
            ball_patches = []                               # List of ball patches.
            position_container = self._container._pos
            r_container = self._r_container
            outline_container = plt.Circle(position_container, r_container,
                                        ec="black", fill=False)

            if random_color:
                for i, ball in np.ndenumerate(self._balls):     # Create tuple pairs.
                    position_ball, r_ball = ball._pos, ball._radius
                    if i != 0:  ball_patches.append(plt.Circle(position_ball, 
                                            r_ball, ec="black", fc=tuple((np.random.rand(),
                                            np.random.rand(), np.random.rand()))))
                    else:       ball_patches.append(plt.Circle(position_ball, r_ball,
                                            ec="black", fc="yellow"))
            
            self._ball_patches = ball_patches
            self._outline_container = outline_container
            
        if draw_current_state:
            init_patches(); plt.figure(num="Current Simulation State")
            axes = plt.axes(xlim=(-self._r_container, self._r_container),
                            ylim=(-self._r_container,self._r_container),
                                aspect="1", animated=True)
            axes.add_patch(self._container_outline)
            for patch in self._ball_patches:    axes.add_patch(patch)
            plt.show()
            
        for i in range(0, self._N_balls):
            self._ball_patches[i].center = self._balls[i].pos()

    def brownian_investigation(self, init=False, radius=5, mass=10, tracer=True):
        """ brownian_investigation | Brownian motion investigation, with ball "0".
            Position of ball 0 is initialised to (0,0). Remaining balls are
            randomly distributed in the container.
                < PARAMETERS >
                -> init(boolean, optional): Creates the investigation.
                -> radius(float, optional): Radius of ball being investigated.
                -> mass(float, optional): Mass of ball being investigated.
                -> tracer(boolean, optional): Draws out Brownian path travelled
                                              in animation for Ball 0.
                RETURNS
                    (plt.Line2D): path travelled between collisions.
        """
        if init:
            self._balls[0].set_pos(np.array([0.0,0.0]))
            self._balls[0].set_radius(radius)
            self._balls[0].set_mass(mass)
            self.gen_random_positions(start=1)
        
        if tracer:
            trace = plt.Line2D(xdata=[self._balls[0]._pos[0], self._brownian[-1][0]],
                            ydata=[self._balls[0]._pos[1], self._brownian[-1][1]],
                            color="black", alpha=0.6,lw=0.7)
            return trace

    def record_property(self, pressure=False, pressure_moyen=False,
                        temperature=False, temperature_moyen=False,
                        speed=False, KE=False, distance_absolute=False,
                        distance_relative=False, brownian=False, df=False,
                        dataset=False):
        """ record_property | Records various properties of the system.
        PROPERTIES
        -   pressure | Record the pressure for every _(*)_ container collisions.
            RAISES
                Exception:
                    IndexError: Insufficient number of collisions.
            RECORDS
            pd.DataFrame [pressure, t]
                CONTAINS
                -> pressure(float): Systemic pressure.
                -> t(float): Time.
                
        -   pressure_moyen | Records the average systemic pressure.
                N.B: Must enable Simulation.run(pressure = True).
                RAISES
                    IndexError: Insufficient number of collisions.
                RECORDS
                    (float): average steady-state systemic pressure.
                    
        -   temperature | Records systemic temperature, ∀ collisions.
                RECORDS
                -> pd.DataFrame [T, t, collision]
                    CONTAINS
                    -> T(float): Systemic temperature.
                    -> t(float): Time.
                    -> collision(int): Collision number.
                    
        -   temperature_moyen | Records the average temperature of the system.
                RECORDS
                    (float): Average systemic temperature.
                    
        -   speed | Record speed for all balls, for all collisions.
                RECORDS
                    (list(float)): Speed of balls, for all collisions.
            
        -   KE | Record the total_KE of the system, for all collisions.
                RECORDS
                -> pd.DataFrame [KE, t, collision]
                    CONTAINS
                    -> KE(float): Systemic kinetic energy.
                    -> t(float): Time.
                    -> collision(int): Collision number.
                    
        -   distance_absolute | Writes distances from origin of all balls
                                        for all collisions.
                RECORDS
                        (list(float)): Ball distances from the centre of container
                        for all collisions.

        -   distance_relative | Writes the relative distances between all
                                        balls for all collisions.
                RECORDS
                    (list(float)): Relative distances between all pairs of balls 
                                   for all collisions.
                    
        -   brownian | Records data prerequisite for Brownian motion
                                investigation.
                RECORDS
                -> pd.DataFrame [x, y, t, collision]
                    CONTAINS
                    -> x(float): x co-ordinate of ball.
                    -> y(float): y co-ordinate of ball.
                    -> t(float): Collision time.
                    -> collision(float): Collision number.
                    
        -   dataset | Write simulation data into pd.DataFrame.
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
        if pressure:    # Revisit and alter.
            self._pressure, N_collisions = [], 50   # Define (*)'s value.
            if not isinstance(self._dp_container, np.ndarray):
                self._dp_container = np.array(self._dp_container)

            try: max_t = self._dp_container[-1,1]   # Error handling.
            except IndexError:                      # No container collisions.
                print("Insufficient number of collisions for pressure calculation.")
                self._pressure = np.nan             # Undefined pressure.
                return

            min_t, start = max_t * 0.2, 0              # Select portion of data.
            while self._dp_container[start, 1] <= min_t: # Find new start index.
                start += 1
                if start == len(self._dp_container - 1):
                    print("Insufficient number of collisions.")
                    self._pressure = np.nan
                    return

            start += (len(self._dp_container) - start) % N_collisions
            dp_new = self._dp_container[start:,:]
            N_pressure = int(len(dp_new)/N_collisions)

            for i in range(N_pressure):
                index = i * N_collisions         # <Pressure> = Σdp / T (period)
                pressure = np.sum(dp_new[index:index+N_collisions-1,0] / (dp_new\
                                [index+N_collisions-1,1]-dp_new[index,1])*2*np.pi\
                                *self._r_container)
                time = (dp_new[index+N_collisions-1,1] + dp_new[index,1]) / 2
                self._pressure.append([pressure,time])

            self._pressure = pd.DataFrame(self._pressure,\
                                        columns = ['pressure', 't'])
        
        if pressure_moyen:
            if not isinstance(self._dp_container, np.ndarray):
                self._dp_container = np.array(self._dp_container)

            try: max_t = self._dp_container[-1,1]              # Error handling.
            except IndexError:                        # No container collisions.
                print("Insufficient number of collisions for pressure calculation.")
                self._pressure = np.nan
                return

            min_t = max_t * 0.2
            start = 0
            while self._dp_container[start, 1] <= min_t: # Find new start index.
                start += 1
                if start == len(self._dp_container - 1):
                    print("Insufficient number of collisions.")
                    self._pressure_moyen = np.nan
                    return
            min_t = self._dp_container[start,1]  # <Pressure> = Σdp / T (period)
            self._pressure_moyen = np.sum(self._dp_container[start:,0]) / ((max_t\
                                        - min_t) * 2 * np.pi * self._r_container)

        if temperature:
            KE = np.zeros(self._N_balls)
            for i, ball in np.ndenumerate(self._balls):
                KE[i] = 0.5 * ball._mass * np.linalg.norm(ball._vel) ** 2
            temperature = (np.sum(KE)) / (self._N_balls * spc.Boltzmann)

            self._temperature.append([temperature, self._utc, self._N_collisions])

            if len(self._temperature) == self._collisions + 1:  # Record data.
                self._temperature = pd.DataFrame(self._temperature,
                                                columns=["T", "t", "collision"])
                
        if temperature_moyen: 
            self._temperature_moyen = np.mean(self._temperature["T"])

        if speed:
            for ball in self._balls:
                self._speed.append(np.sqrt(np.dot(ball._vel,ball._vel)))

        if np.any(KE):      # Have no idea, but wants np.any used.
            total_KE = np.sum([0.5 * self._m_balls * np.dot(ball._vel,ball._vel)\
                            for ball in self._balls])
            self._KE.append([total_KE, self._utc, self._N_collisions])

            if len(self._KE) == self._collisions + 1:
                self._KE = pd.DataFrame(self._KE, columns=['KE','t','collision'])

        if distance_absolute:
            for ball in self._balls:
                self._distance_absolute.append(np.sqrt(np.dot(ball._pos,
                                                              ball._pos)))

        if distance_relative:
            for _, pair in enumerate(self._pairs):      # Consider changing to np.
                ball_A, ball_B = pair[0], pair[1]
                distance_relative = np.sqrt(np.dot(self._balls[ball_A]._pos - \
                    self._balls[ball_B]._pos, self._balls[ball_A]._pos - \
                    self._balls[ball_B]._pos))
                self._distance_relative.append(distance_relative)         

        if brownian:
            if df == False:
                self._brownian.append(np.array([self._balls[0]._pos[0],\
                    self._balls[0]._pos[1], self._utc, self._N_collisions]))
            else:       # Once data completes, add to pd.DataFrame
                self._brownian = pd.DataFrame(self._brownian,\
                                columns = ['x','y','t','collision'])
                
        if dataset:
            for i, ball in enumerate (self._balls):
                j = self._N_collisions * self._N_balls
                self._dataset[i+j, 0] = i
                self._dataset[i+j, 1] = ball._mass
                self._dataset[i+j, 2] = ball._pos[0]
                self._dataset[i+j, 3] = ball._pos[1]
                self._dataset[i+j, 4] = ball._vel[0]
                self._dataset[i+j, 5] = ball._vel[1]
                self._dataset[i+j, 6] = self._N_collisions
                self._dataset[i+j, 7] = self._utc

                if (np.sqrt(np.dot(ball._pos,ball._pos)) - self._r_container + \
                    ball._radius <= 10e-10):    # Takes into account floating point.
                    self._dataset[i+j, 8] = True
                else: self._dataset[i+j, 8] = False

            if self._N_collisions == self._collisions:
                self._dataset = pd.DataFrame(self._dataset, columns=['ball','mass',\
                                    'x','y','v_x','v_y','collision','t','container'])

    def record_data_states(self, distance_absolute = False,
                           distance_relative = False, speed = False, KE = False,
                           test_temperature = False, temperature = False,
                           dataset = False, pressure=False, test_pressure=False):
        """ record_data_states | Records simulation datasets.
            < PARAMETERS >
                -> distance_absolute(boolean, optional): If True, writes dataset
                    for all distances of ball to O, for all collisions.
                -> distance_relative(boolean, optional): If True, writes dataset
                    for all relative distances between balls, for all collisions.
                -> speed(boolean, optional): If True, does the same as above for
                    speeds of all balls, for all collisions.
                -> KE(boolean, optional): If True, same as above, for system KE.
                -> temperature(boolean, optional): If True, same as above, for
                    the average temperature of the system (not ∀ collisions).
                -> dataset(boolean, optional): If True, same as above, for all
                    simulation information.
                -> pressure(boolean, optional): If True, writes dataset for the
                                                average systemic pressure.
                -> test_pressure(boolean, optional): If True, writes dataset of
                            pressure for every _(*)_ container collisions.
        """
        if distance_absolute:                   self.record_property(distance_absolute=True)
        if distance_relative:                   self.record_property(distance_relative=True)
        if speed:                               self.record_property(speed=True)
        if KE:                                  self.record_property(KE=True)
        if temperature or test_temperature:     self.record_property(temperature=True)
        if dataset:                             self.record_property(dataset=True)
        if pressure:                            self.record_property(pressure_moyen=True)
        if test_pressure:                       self.record_property(pressure=True)

    def run(self,
        collisions = 10,            # The number of collisions.
        time = 0.001,               # Time period between collisions.

        pressure = False,           # Enable pressure calculations.
        temperature = True,         # Enable temperature calculations.
        KE = False,                 # Enable kinetic energy calculations.
        speed = False,              # Enable speed calculations.
        brownian = False,           # Enable Brownian motion inveestigation.
        dataset = False,            # Enable complete dataset.
        test_pressure = False,      # Enable test pressure calculations.
        test_temperature = False,   # Enable test temperature calculations.
        distance_absolute = False,  # Use ball distance from the origin.
        distance_relative = False,  # Use relative distances between balls.
        progress_bar = True,            # Enable progress bar animation in terminal.
        animate = False,            # Enables animation.
        random_color = True,        # Assigns random colors to balls.
        sim_title = "1. Simulation" # Give the simulation title.
    ):
        """ run | Runs the 2D simulation of colliding particles within the container.
        PARAMETERS
        -> collisions (int, optional): number of collisions in the simulation.
        -> time (float, optional): time period (s) between animation frames.

        -> pressure (boolean, optional): records pressure for every __ collisions
                                         with the wall of the container.
        -> temperature (boolean, optional): records temperature.
        -> KE (boolean, optional): records the system KE for every collision.
        -> speed (boolean, optional): records speed of all balls in all collisions.
        -> brownian (boolean, optional): records data for Brownian motion.
        -> dataset (boolean, optional): records dataset of simulation information.
        -> distance_absolute (boolean, optional): records the distances of all
                                    balls from the origin, in all collisions.
        -> distance_relative (boolean, optional): records the relative distances
                                    between all the balls, in all collisions.
        -> progress_bar (boolean, optional): displays a progress bar.
        -> sim_title (str, optional): Provide the title for the plot.

        RETURNS
            gloss(dict): Glossary of needed datasets.
        """
        self._utc = 0
        self._queue = hd.heapdict()        # Create priority queue.
        self._collisions = collisions
        self._animate = animate
        self._has_brownian = brownian
        self._time = time

        print("Starting {self._N_balls} balls, r_balls = {self._r_balls}, \
            speed range = {self._speed_range}, r_container = {self._r_container},\
            {self._collisions} collisions.")
        
        if animate:                             # Initialise animation.
            self.draw(init=True)
            plt.figure(num="Thermodynamic Simulation, Animated")
            plt.rcParams.update(plt.rcParamsDefault)
            ax = plt.axes(xlim=(-self._r_container,self._r_container),
                          ylim=(-self._r_container,self._r_container),aspect="1",
                          xticks=(),yticks=())
            ax.add_patch(self._outline_container)
            ax.set_title("Simulation Animation", **titleFont)
            ax.set_xlabel("displacement in x / m", **axesFont)
            ax.set_ylabel("displacement in y / m", **axesFont)
            # ax.set_xticks(np.linspace(-10,10,3),**axesFont)
            # ax.set_yticks(np.linspace(-10,10,3),**axesFont)
            self._axes = ax
            for ball_patch in self._ball_patches: ax.add_patch(ball_patch)
            plt.pause(time)

        if dataset:         self._dataset = np.zeros((self._N_balls *
                                                    (self._collisions + 1), 9))
        if brownian:        self.record_property(brownian=True)
        if progress_bar:    self._time_epoch = tm.time()

        self.record_data_states(temperature=temperature, speed=speed, KE=KE,
                                distance_absolute=distance_absolute,
                                distance_relative=distance_relative,
                                dataset=dataset)
        self.next_collision(init=True)          # Run the first collision.
        self.record_data_states(temperature=temperature, speed=speed, KE=KE,
                                distance_absolute=distance_absolute,
                                distance_relative=distance_relative,
                                dataset=dataset)

        for i in range(2, collisions + 1):      # For subsequent collisions:
            if progress_bar:    progress(self._time_epoch, i, collisions)
            self.next_collision(init=False)     # Run subsequent collisions.
            self.record_data_states(temperature=temperature, speed=speed, KE=KE,
                                    distance_absolute=distance_absolute,
                                    distance_relative=distance_relative,
                                    dataset=dataset)
            # plt.savefig("animation/"+str(i)+".png",dpi=1000)

        if animate:         plt.show()
        if temperature:     self.record_property(temperature_moyen=True)
        self.record_data_states(pressure=pressure, test_pressure=test_pressure)
        if brownian:    self.record_property(brownian=True, df=True)

        gloss = self.glossary(temperature=temperature,
                              test_temperature=test_temperature,
                              pressure=pressure, test_pressure=test_pressure,
                              speed=speed, KE=KE,
                              distance_absolute=distance_absolute,
                              distance_relative=distance_relative,
                              dataset=dataset, brownian=brownian)

        print("Ending {self._N_balls} balls, r_balls = {self._r_balls},\
            speed range = {self._speed_range}, r_container = {self._r_container},\
            {self._collisions} collisions.")
        return gloss

def progress(start_time, iterations, all_iterations, description="Collision"):
    """ progress | A progress bar to show the progression of an operation.
        < PARAMETERS >
        -> start_time: The start time of the collisions.
        -> iterations: The number of iterations performed.
        -> all_iterations: The total number of iterations:
        -> description: Description.
    """
    current_time = tm.time()
    n_steps = 50
    t_elapsed = current_time - start_time
    t_remaining = t_elapsed * all_iterations / iterations - t_elapsed
    t_elapsed_str = tm.strftime("%H:%M:%S", tm.gmtime(t_elapsed))
    t_remaining_str = tm.strftime("%H:%M:%S", tm.gmtime(t_remaining))

    percentage = round(iterations / all_iterations * 100)
    n_blocks = int(np.floor(percentage / (100/n_steps)))
    blocks = n_blocks * "\u2588" + (n_steps - n_blocks) * " " # Alter
    description_s = f"{description}:"

    if t_elapsed == 0:  it_s = "0.00"   # Calculating iterations per second.
    else:               it_s = round(iterations / t_elapsed, 2)

    progress = f"\r[{blocks}] {percentage}% complete. | {description_s} {iterations}/{all_iterations} | Time: {t_elapsed_str}/{t_remaining_str} | Speed: " + str(np.round(it_s, 1)) + " collisions.s¯¹ ||"

    if iterations == all_iterations:    progress += "\n"
    print ("\033[A \033[A");            sys.stdout.write(progress)

class Event(tuple):
    """ EVENT CLASS |
        A tuple of 5 elements (ball_A, ball_B, count_A, count_B, dt).
        < PARAMETERS >
        -> ball_A (int): The first ball in impending collision.
        -> ball_B (int): The second ball in impending collision.
        -> count_A (int): The number of collisions the first ball
                encountered prior to this impending collision calculation.
        -> count_B (int): The number of collisions the second ball
                encountered prior to this impending collision calculation.
        -> dt (float): The global time this collision will happen on.
        RETURNS
            pair (list(int)): a list of the two balls or a ball and container
                              involved in the collision.
            ball_A (int): the index of first ball in impending collision.
            ball_B (int): the index of second ball in impending collision.
            count_A (int): the number of collisions the first ball encountered
                           prior to this impending collision calculation.
            count_B (int): the number of collisions the second ball encountered
                           prior to this impending collision calculation.
            dt (float): the global time(step) this collision will happen on.
    """
    def pair(self):         return [self[0], self[1]]   # Returns pair object.
    def ball_A(self):       return self[0]              # Returns Ball A.
    def ball_B(self):       return self[1]              # Returns Ball B.
    def count_A(self):      return self[2]              # Returns Ball A count.
    def count_B(self):      return self[3]              # Returns Ball B count.
    def dt(self):           return self[4]              # Returns minimum time.
    
class Plot:
    def plot_simulation_animation(m_balls_animation=5e-26, r_balls_animation=0.2,
                                N_balls_animation=100, collisions_animation=500,
                                r_container_animation=10,
                                random_speed_range_animation=500):
        """ plot_simulation_animation | Plot the animation of the simulation.
                < PARAMETERS >
                -> m_balls_animation(float, optional): Ball mass in animation.
                -> r_balls_animation(float, optional): Ball radius in animation.
                -> N_balls_animation(int, optional): Number of animation balls.
                -> collisions_animation(int, optional): Number of collisions.
                -> r_container_animation(float, optional): Radius of container.
                -> random_speed_range_animation(int, optional): Range of speeds.
        """
        print("(1) Commencing the animated Simulation...")
        simulation_animation = Simulation(m_balls=m_balls_animation, 
            r_balls=r_balls_animation, N_balls=N_balls_animation, r_container=
            r_container_animation, random_speed_range=random_speed_range_animation)
        simulation_animation.run(collisions=collisions_animation, animate=True)

    def generate_dataset(folder_name="data", m_balls=5e-26, N_balls=[10,20],
                         r_balls=[2,4], r_containers=[100,200], collisions=50,
                         speed_ranges=[500,1000]):
        """ generate_dataset | Generates a dataset
            < PARAMETERS >
            -> folder_name(str): Designated name of folder holding datasets.
            -> m_balls(float, optional): Mass of the dataset balls.
            -> N_balls(list(int), optional): List of number of balls to create.
            -> r_balls(list(float), optional): List of ball radii to create.
            -> r_containers(list(float), optional): List of container radii.
            -> collisions(int, optional): Number of collisions to occur.
            -> speed_ranges(list(float), optional): List of ball speed ranges.
        """
        print("(2) Generating exportable CSV datasets...")
        
        DATA_PATH = os.path.join(os.getcwd(), folder_name)       # Create path.
        if not os.path.exists(DATA_PATH): os.makedirs(DATA_PATH) # Make folder.
        
        for r_ball in r_balls:
            for N_ball in N_balls:
                for r_container in r_containers:
                    for random_speed_range in speed_ranges:
                        fname = f"dataset_{N_ball}_{r_ball}_{r_container}_\
                                 {m_balls}_{random_speed_range}_{collisions}.csv"
                        FILE_PATH = os.path.join(DATA_PATH, fname)
                        if os.path.exists(FILE_PATH):
                            print(f"Existing Dataset: {N_ball} balls,\
                                    r_ball = {r_ball}, r_container = {r_container},\
                                    max speed = {random_speed_range},\
                                    {collisions} collisions")
                            continue
                        
                        else:
                            sim_gen_dataset = Simulation(m_balls=m_balls, r_balls=r_ball,
                                         r_container=r_container, N_balls=N_ball,
                                         random_speed_range=random_speed_range)
                            dataset = sim_gen_dataset.run(collisions=collisions, 
                                                          dataset=True)["dataset"]
                            dataset.to_csv(FILE_PATH)
                            print(f"Generated Dataset: {FILE_PATH}")
                            
    def plot_distance(init=True, absolute=True, relative=True):
        """ plot_distance | Plot relative and absolute distance histograms.
                < PARAMETERS >
                -> init(boolean, optional): Check if this is the initialisation.
                -> absolute(boolean, optional): Plot absolute distances.
                -> relative(boolean, optional): Plot relative distances.
        """
        print("(3) Commencing the distance histograms...")
        simulation_distance = Simulation(N_balls=50, r_container=10,r_balls=0.2,
                                m_balls=5e-26, random_speed_range=500)
        sim = simulation_distance.run(collisions=5000, KE=True,
                                        distance_absolute=True,
                                        distance_relative=True)
        distance_absolute_plot = sim["distance from centre"]
        distance_relative_plot = sim["relative distance"]
        
        if absolute:
            # plt.figure(num="Absolute Distance from Origin")
            sns.set(context="paper", style="darkgrid", palette="muted")
            sns.displot(distance_absolute_plot, kde=False, kind="hist", **histStyle)
            plt.title("Distribution of Ball Absolute Distances", **titleFont)
            plt.xlabel("Ball Distance from Origin $/m$", **axesFont)
            plt.ylabel("Probability Density $/m^{-1}$", **axesFont)
            plt.xticks(**ticksFont); plt.yticks(**ticksFont) 
            plt.tight_layout(); plt.show()
            
        if relative:
            # plt.figure(num="Relative Distance of Ball Pairs")
            sns.set(context="paper", style="darkgrid", palette="muted")
            sns.displot(distance_relative_plot, kind="hist", **histStyle)
            plt.title("Distribution of Ball Relative Distances", **titleFont)
            plt.xlabel(r"Distance between Balls $/m$", **axesFont)
            plt.ylabel(r"Probability Density $/m^{-1}$", **axesFont)
            plt.xticks(**ticksFont); plt.yticks(**ticksFont) 
            plt.tight_layout(); plt.show()

    def plot_conservation(init=True, pressure=False, KE=False, p=False):
        """ plot_conservation | Plots conservation graphs over time.
                < PARAMETERS >
                -> init(boolean, optional): Initialisation method.
                -> pressure(boolean, optional): Enable pressure graph.
                -> KE(boolean, optional): Enable KE graph.
                -> p(boolean, optional): Enable momentum graph.
        """
        print("(4) Commencing the conservation law plots...")
        if init:
            print("Initialisation")
        if pressure:
            print("Pressure Conservation...")
            print("Starting Simulation for Steady State Pressure Distribution")
            
            m_ball_P = 5e-26
            N_ball_P = 10
            r_ball_P = 0.2
            r_container_P = 1
            random_speed_range_P = 500
            collisions_P = 50000

            sim_test_P = Simulation(
                N_balls=N_ball_P,
                r_container=r_container_P,
                r_balls=r_ball_P,
                m_balls=m_ball_P,
                random_speed_range=random_speed_range_P,
            )
            param_test_P = sim_test_P.run(collisions=collisions_P, test_pressure=True)

            pressure_test = param_test_P["pressure"]

            # Fitting the Steady State Pressure to a Normal Distribution
            fit_norm = sp.stats.norm.fit(pressure_test["pressure"])
            ks_norm = sp.stats.kstest(pressure_test["pressure"], "norm", fit_norm)
            pdf_norm = sp.stats.norm(*fit_norm)

            print(f"Normal Distribution: KS Statistic = \
                {round(ks_norm[0],6)}, p = {round(ks_norm[1],3)}")


            arr_norm = np.linspace(min(pressure_test["pressure"]),\
                max(pressure_test["pressure"]), 1000)

            legend_norm = r"$\mu$ = %s Pa, $\sigma$ =  %s Pa" % (
                float("%.2g" % fit_norm[0]),
                float("%.1g" % fit_norm[1]))

            print("Plotting Graph 1 of 4")

            plt.figure(num="Steady State Pressure Distribution")
            sns.set(context="paper", style="darkgrid", palette="muted")
            sns.distplot(pressure_test["pressure"], kde=False, norm_hist=True, rug=True)
            plt.plot(arr_norm, pdf_norm.pdf(arr_norm), label=legend_norm, lw=2)
            plt.title("Steady State Pressure Distribution")
            plt.xlabel(r"Pressure $/Pa$")
            plt.ylabel(r"Probability Density $/Pa ^{-1}$")
            plt.legend()
            plt.tight_layout()
            plt.show()

        if KE:
            print("Conservation of Energy...")
        if p:
            print("Conservation of Pressure...")
    
    def plot_ideal_gas_law(m_ball=5e-26, r_ball=0.1, collisions=500,
                           r_containers=[50,60,70,80,90,100,110,120,130,140,150],
                           N_balls=[50,100], random_speed_ranges=[500,1000]):
        """ plot_ideal_gas_law | Plots a graph of Pressure against Temperature
                                 according to the ideal gas law.
        """
        def ideal_gas_law(V, N, T):
            """ ideal_gas_law | Use ideal gas law, returning P.
            < PARAMETERS >
            -> V (float): Volume of the container.
            -> N (int): Number of gas particles.
            -> T (float): The temperature of the system.
            RETURNS
                (float): The pressure of the system.
            """
            return V ** -1 * N * spc.Boltzmann * T

        def run_ideal_gas_simulation(params):
            """ run_ideal_gas_simulation | Run 2D simulation of Ball collisions.
                    < PARAMETERS >
                    -> parameters(list[N_ball, r_container, random_speed_range,
                                  speed]): 
                        N_ball(int): The number of Ball objects.
                        r_container(float): The radius of Container.
                        random_speed_range(float): The speed range generated.
                        speed(list(np.ndarray[float])): Ball velocity list.
            """
            N_ball, r_container, random_speed_range, speed = parameters[0], 
            parameters[1], parameters=[2], parameters[3], parameters[4]
            
            simulation_ideal_gas = Simulation(N_balls=N_ball, r_balls=r_ball,
                                              r_container=r_container, 
                                              random_speed_range=random_speed_range,
                                              m_balls=m_ball)
            simulation_ideal_gas.set_vel_ball(speed)
            sim = simulation_ideal_gas.run(collisions=collisions, pressure=True,
                                           temperature=True, progress=True)
            return sim
        
        parameter_list, pressures, temperatures, speeds = [], [], [], []
        volumes = [np.pi*r**2 for r in r_containers]
        
        if __name__ == "__main__":
            print("Generating speeds...")
            for N in N_balls:
                for random_speed_range in random_speed_ranges:
                    speeds.append(Simulation.gen_random_velocities(N_balls=N, 
                                                        max_speed=random_speed_range))
            i = 0
            for N in N_balls:
                for random_speed_range in random_speed_ranges:
                    for r_container in r_containers:
                        parameter_list.append([N, random_speed_range, r_container,
                                               speeds[i]])
                    index += 1
            print("Commencing simulations...")
            start_time = tm.perf_counter()
            with ft.ProcessPoolExecutor() as executor:
                results = executor.map(run_ideal_gas_simulation, parameter_list)
            end_time = tm.perf_counter()
            print(f"Time elapsed: {np.round(end_time-start_time,2)} seconds.")
            
            for i, result in enumerate(results):
                pressures.append(result["average temperature"])
                if i % len(volumes) == 0:
                    temperatures.append(result["average temperature"])
                    
            array_fit = np.linspace(volumes[0], volumes[-1], 1000)
            print("Plotting graph..."); j = 0
            
            plt.figure(num="Ideal Gas Law over Varying Temperatures")
            sns.set(context="paper", style="darkgrid", palette="muted")
            
            for N in N_balls:
                for random_speed_range in random_speed_ranges:
                    legend = f"N = {N}, T = %s K" % (np.float("%.3g" 
                                                     % temperatures[j]))
                    plt.plot(array_fit, ideal_gas_law(array_fit, N,
                                        temperatures[j]), label=legend)
                    j += 1
            
            for i, _ in enumerate(temperatures):
                pressures_temp = []
                for j, _ in enumerate(volumes):
                    pressures_temp.append(pressures[i*len(volumes) + j])
                plt.plot(volumes, pressures_temp, "o", mec="white", mew=0.5)
            
            plt.title("Ideal Gas Law")
            plt.xlabel("Volume /$m^2$")
            plt.ylabel("Pressure /Pa")
            plt.legend(); plt.tight_layout(); plt.show()

    def plot_maxwell_boltzmann(m_ball=5e-26, r_container=5, N_balls=300,
                               r_ball=0.1, collisions = 10000, random_speed_range=500,
                               range_N_balls = np.linspace(100,300,101), save=False,
                               range_N_collisions = np.linspace(125, 10000, 80),
                               iter_balls=True,iter_collisions=True):
        """ plot_maxwell_boltzmann | Plots a Maxwell-Boltzmann distribution.
        """
        range_N_balls = (int(N) for N in range_N_balls)
        range_N_collisions = (int(N) for N in range_N_collisions)
        temperature_list = []
        def distribution_maxwell_boltzmann(v, m, T):
            """ distribution_maxwell_boltzmann | Returns value of distribution.
                    < PARAMETERS >
                    -> v(np.ndarray(float)): Array of gas particle speeds. 
            """
            return (m*v)/(spc.Boltzmann*T)*np.exp((-0.5*m*v**2)/(spc.Boltzmann*T))
        
        if iter_balls:
            for N_ball in range_N_balls:
                simulation_maxwell_boltzmann = Simulation(N_balls=N_ball,
                            r_container=r_container, m_balls=m_ball, r_balls=r_ball,
                            random_speed_range=random_speed_range)
                sim = simulation_maxwell_boltzmann.run(collisions=collisions, speed=True,
                                                    temperature=True)
                speeds = np.array(sim["speed"])
                temperature = sim["average temperature"]
                
                array_maxwell_boltzmann = np.linspace(np.amin(speeds),np.amax(speeds),10000)
                
                n, bins, patches = plt.hist(x = speeds, bins = 30, **histStyle)
                v = (bins[:-1] + bins[1:]) / 2
                curve_maxwell_boltzmann = opt.curve_fit(distribution_maxwell_boltzmann, n, v)
                
                # curve_maxwell_boltzmann = distribution_maxwell_boltzmann(
                    #array_maxwell_boltzmann, temperature, m_ball)
                
                font = fnt.FontProperties(family='C059', weight='bold',
                                          style='normal', size=8)

                x = np.linspace(bins[0], bins[-1], 1000)
                params, cov_params = opt.curve_fit(distribution_maxwell_boltzmann, v, n)
                
                sns.set(context="paper", style="darkgrid", palette="muted")
                # sns.displot(speeds, label="Simulation Data", bins=30, kde=False)
                
                plt.plot(x, distribution_maxwell_boltzmann(x, *params), 
                        label="Maxwell-Boltzmann Distribution", lw=2)
                
                plt.title("2D Maxwell-Boltzmann Distribution", **titleFont)
                plt.xlabel(r"Speed /$m s^{-1}$ ", **axesFont)
                plt.ylabel(r"Probability Density /$m^{-1} s$", **axesFont)
                plt.xticks(**ticksFont); plt.yticks(**ticksFont)
                plt.legend(prop=font)
                plt.tight_layout()
                plt.show()

                temperature_list.append(temperature)
                if save:
                    np.savetxt(f"data/speeds_m_{m_ball}_r_c_{r_container}_N_{N_ball}_r_b_{r_ball}_N_c_{collisions}_r_s_r_{random_speed_range}.txt", speeds, fmt="%s")
                    np.savetxt(f"data/arr_m_{m_ball}_r_c_{r_container}_N_{N_ball}_r_b_{r_ball}_N_c_{collisions}_r_s_r_{random_speed_range}.txt", array_maxwell_boltzmann, fmt="%s")
                    np.savetxt(f"data/MBD_m_{m_ball}_r_c_{r_container}_N_{N_ball}_r_b_{r_ball}_N_c_{collisions}_r_s_r_{random_speed_range}.txt", curve_maxwell_boltzmann, fmt="%s")
                    np.savetxt(f"data/temperature_m_{m_ball}_r_c_{r_container}_N_{N_ball}_r_b_{r_ball}_N_c_{collisions}_r_s_r_{random_speed_range}.txt", temperature_list, fmt="%s")

        if iter_collisions:
            temperature_list = []
            for collisions in range_N_collisions:
                simulation_maxwell_boltzmann = Simulation(N_balls=N_ball,
                            r_container=r_container, m_balls=m_ball, r_balls=r_ball,
                            random_speed_range=random_speed_range)
                sim = simulation_maxwell_boltzmann.run(collisions=collisions, speed=True,
                                                    temperature=True)
                speeds = np.array(sim["speed"])
                temperature = sim["average temperature"]
                
                array_maxwell_boltzmann = np.linspace(np.amin(speeds),np.amax(speeds),
                                                    10000)
                curve_maxwell_boltzmann = distribution_maxwell_boltzmann(
                    array_maxwell_boltzmann, temperature, m_ball)

                font = fnt.FontProperties(family='C059', weight='bold',
                                          style='normal', size=8)

                sns.set(context="paper", style="darkgrid", palette="muted")

                sns.displot(speeds, label="Simulation Data", bins=30, kde=False)
                plt.plot(array_maxwell_boltzmann, curve_maxwell_boltzmann, 
                        label="Maxwell-Boltzmann Distribution", lw=2)
                plt.title("2D Maxwell-Boltzmann Distribution", **titleFont)
                plt.xlabel(r"Speed /$m s^{-1}$ ", **axesFont)
                plt.ylabel(r"Probability Density /$m^{-1} s$", **axesFont)
                plt.xticks(**ticksFont); plt.yticks(**ticksFont)
                plt.legend(prop=font)
                plt.tight_layout()
                plt.show()


                temperature_list.append(temperature)
                if save:
                    np.savetxt(f"data/c/speeds_m_{m_ball}_r_c_{r_container}_N_{N_ball}_r_b_{r_ball}_N_c_{collisions}_r_s_r_{random_speed_range}.txt", speeds, fmt="%s")
                    np.savetxt(f"data/c/arr_m_{m_ball}_r_c_{r_container}_N_{N_ball}_r_b_{r_ball}_N_c_{collisions}_r_s_r_{random_speed_range}.txt", array_maxwell_boltzmann, fmt="%s")
                    np.savetxt(f"data/c/MBD_m_{m_ball}_r_c_{r_container}_N_{N_ball}_r_b_{r_ball}_N_c_{collisions}_r_s_r_{random_speed_range}.txt", curve_maxwell_boltzmann, fmt="%s")
                    np.savetxt(f"data/c/temperature_m_{m_ball}_r_c_{r_container}_N_{N_ball}_r_b_{r_ball}_N_c_{collisions}_r_s_r_{random_speed_range}.txt", temperature_list, fmt="%s")