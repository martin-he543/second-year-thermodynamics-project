"""
    SIMULATION MODULE (import simulation as sm) | Created by 𝑀𝒶𝓇𝓉𝒾𝓃 𝒜. 𝐻𝑒, 2022.11.20
    For the simulation of elastic collisions of balls with others, and container.
"""
    #%% Simulation Initialisation ✔️
import ball as bl
import event as ev
import numpy as np
import scipy.constants as spc
import matplotlib.pyplot as plt
import pandas as pd
import itertools as it
import heapdict as hd
import time as tm
import sys, os, random

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
        self._ball = [] #!!!!!!! EDIT

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
        self._N_collisions = 0
                
        for _ in range(0, N_balls):
            self._ball.append(bl.Ball(radius=r_balls, mass=m_balls))
        if random_position: self.generator_random_position()
        self.generator_random_vel(max_speed=random_speed_range)
        
    #%% Simulation Information ✔️
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
    
    #%% Simulation Property ✔️
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
    speed | Provides the speed for all balls, in all collisions.
        N.B: Must enable Simulation.run(speed = True).
        RETURNS
            speed(list(float)): ball speeds, in all collisions.

    KE | Provides the total system KE, for all collisions.
        N.B: Must enable Simulation.run(KE = True).
        RETURNS
            KE(list(float)): Total systemic KE, at all collisions.

    distance_absolute | Provides distances from the origin for all balls, in all collisions.
        N.B: Must enable Simulation.run(distance_absolute = True).
        RETURNS
            distance_absolute(list(float)): ball distances from the origin, for all collisions.

    distance_relative | Provides relative distance between all balls, in all collisions.
        N.B: Must enable Simulation.run(distance_relative = True).
        RETURNS
            relative(list(float)): ball distances between all possible pairs of balls, in all collisions.

    temperature | Provides the system temperature at all collision times.
        N.B: Must enable Simulation.run(temperature = True).
        RETURNS
            temperature(list(float)): Systemic temperature, at all collision times.
        
    temperature_moyen | Provides the average temperature of the system.
        N.B: Must enable Simulation.run(temperature = True).
        RETURNS
            temperature_moyen(float): the system's average temperature.
                 
    pressure | Provides the pressure values for every __ container collisions.
        N.B: Must enable Simulation.run(pressure = True).
        RETURNS
            pressure(float): system pressure.

    pressure_moyen | Provides the average pressure of the system.
        N.B: Must enable Simulation.run(pressure = True).
        RETURNS
            pressure_moyen(float): the system's average steady-state pressure.

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
    def ball(self):                 return self._ball
    def container(self):            return self._container
    def speed(self):                return self._speed
    def KE(self):                   return self._KE
    def distance_absolute(self):    return self._distance_absolute
    def distance_relative(self):    return self._distance_relative
    def temperature(self):          return self._temperature
    def temperature_moyen(self):    return self._temperature_moyen
    def pressure(self):             return self._pressure
    def pressure_moyen(self):       return self._pressure_moyen
    def brownian(self):             return self._brownian
    def dataset(self):              return self._dataset
    
    def pair_combn(self, container=False):
        if not container:       list_number = list(range(self._N_balls))
        else:                   list_number = list(range(self._N_balls + 1))
        return                  list(it.combinations(list_number, 2))   
  

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
        """
        Initialise next collision time calculations for the first timestep.
        Calculate all possible ball pairs and their respective impending
        collision time.
        Collision times are recorded as an event.Event object.
        All collision events are added into a priority queue for high 
        efficiency selection of next event.
        The priority queue is a binary heap implemented using heapdict.heapdict
        The root node of this priority queue will always be the next immediate 
        event (it has the smallest time value).
        """
        # Calculating all collisions between balls
        for pair in self._pairs:  # All possible ball pair combinations
            ball_A = self._ball[pair[0]]
            ball_B = self._ball[pair[1]]
            dt = ball_A.time_to_collision(ball_B)
            if dt != np.inf:  # Selecting only valid solutions
                self._pq[
                    ev.Event((pair[0], pair[1], ball_A._count, ball_B._count, dt))
                ] = dt  # Adding event to priority queue

        # Calculating collisions between balls and container
        for i, ball in enumerate(self._ball):
            dt = ball.time_to_collision(self._container)
            if dt != np.inf:
                self._pq[ev.Event((i, self._N_balls, ball._count, -1, dt))] = dt
# CHECK
    def update_patch(self):
        """
        Updates the positions of ball patches in animation.
        """
        for i in range(0, self._N_balls):
            self._b_patch[i].center = self._ball[i].pos()
# CHECK
    def trace_brownian(self):
        """
        Draws out the path travelled by ball 0 in animation.

        Returns:
            (matplotlib.pyplot.Line2D): The path travelled between previous and 
                current collision.
        """
        path = plt.Line2D(
            xdata=[self._ball[0]._pos[0], self._brownian[-1][0]],
            ydata=[self._ball[0]._pos[1], self._brownian[-1][1]],
            color="black",
            alpha=0.8,
            lw=1,
        )
        return path
# CHECK
    def collision_time(self):
        """
        Calculates next collision times of the balls that underwent collisions.
        """
        collided_ball = set()
        for event in self._events:  # Events of next collisions
            for collided in event.pair():
                collided_ball.add(collided)

        # Adds collision events to priority queue
        for element in collided_ball:
            if element != self._N_balls:
                # Calculating collisions with container
                dt = self._ball[element].time_to_collision(self._container)
                if dt != np.inf:
                    self._pq[
                        ev.Event(
                            (
                                element,
                                self._N_balls,
                                self._ball[element]._count,
                                -1,
                                dt + self._global_time,
                            )
                        )
                    ] = (dt + self._global_time)

                # Calculating collisions with other balls
                for j in range(self._N_balls):
                    if j != element:
                        # Ensure smaller index comes first
                        if j < element:
                            ball_A = self._ball[j]
                            ball_B = self._ball[element]
                            index_A = j
                            index_B = element
                        else:
                            ball_A = self._ball[element]
                            ball_B = self._ball[j]
                            index_A = element
                            index_B = j
                        dt = ball_A.time_to_collision(ball_B)
                        if dt != np.inf:
                            self._pq[
                                ev.Event(
                                    (
                                        index_A,
                                        index_B,
                                        self._ball[index_A]._count,
                                        self._ball[index_B]._count,
                                        dt + self._global_time,
                                    )
                                )
                            ] = (dt + self._global_time)

            # If container underwent collision
            else:
                for j in range(self._N_balls):
                    dt = self._ball[j].time_to_collision(self._container)
                    if dt != np.inf:
                        self._pq[
                            ev.Event(
                                (
                                    j,
                                    self._N_balls,
                                    self._ball[j]._count,
                                    -1,
                                    dt + self._global_time,
                                )
                            )
                        ] = (dt + self._global_time)
# CHECK
    def init_next_event(self):
        """
        Initialising next event selection, taking into account that multiple 
        collisions might occur at the same time.
        """
        self._events = []  # A list of next events

        min_event = self._pq.popitem()[0]  # Picking next event
        self._min_dt = min_event.dt()
        self._events.append(min_event)

        # Checks if multiple collisions happen at the same time
        while len(self._pq) != 0:
            if self._pq.peekitem()[0].dt() == self._min_dt:
                self._events.append(self._pq.popitem()[0])
            else:
                break
# CHECK
    def next_event(self):
        """
        Selecting the next collision event.
        If the collision count of the ball has increased compared to that of 
        the event, it means that the ball has collided with other balls after 
        the event is calculated, invalidating the event. Such events are 
        discarded.
        """
        self._events = []
        min_event = self._pq.popitem()[0]

        # Checks validity of event
        while len(self._pq) != 0:
            min_A = min_event.ball_A()  # Ball numbers
            min_B = min_event.ball_B()
            if min_B == self._N_balls:  # Container collision
                if min_event.count_A() != self._ball[min_A]._count:
                    min_event = self._pq.popitem()[0]  # Picks next event
                else:
                    break
            else:  # Collision with other balls
                if (
                    min_event.count_A() != self._ball[min_A]._count
                    and min_event.count_B() != self._ball[min_B]._count
                ):
                    min_event = self._pq.popitem()[0]  # Picks next event
                else:
                    break
                # check for invalidated collision

        self._min_dt = min_event.dt()
        self._events.append(min_event)

        # Checks if there are other events with the same collision time
        while len(self._pq) != 0:
            next_event = self._pq.peekitem()[0]
            if next_event.dt() == self._min_dt:
                next_A = next_event.ball_A()  # Ball numbers
                next_B = next_event.ball_B()
                if next_B == self._N_balls:  # Container collision
                    if next_event.count_A() == self._ball[next_A]._count:
                        self._events.append(self._pq.popitem()[0])
                    else:
                        break
                else:  # Collision with other balls
                    if (
                        next_event.count_A() == self._ball[next_A]._count
                        and next_event.count_B() == self._ball[next_B]._count
                    ):
                        self._events.append(self._pq.popitem()[0])
                    else:
                        break
            else:
                break

    def move_balls(self):
        """
        Moves balls to the timestep of next collision.
        """
        for ball in self._ball:
            ball.move(self._min_dt - self._global_time)

    def collide_balls(self, pressure, test_pressure, brownian):
        """
        Collides balls, changing their velocities.

        Parameters:
            pressure (boolean): If True, pressure data is recorded.
            test_pressure (boolean): If True, pressure data is recorded.
            brownian (boolean): If True, data for Brownian Motion is recorded.
        """
        record = False
        for event in self._events:
            ball_1 = event.ball_A()
            ball_2 = event.ball_B()
            if ball_2 == self._N_balls:  # Container collision
                self._ball[ball_1].collide(self._container)
                self._N_container_collisions += 1
                if pressure or test_pressure:  # Appends change in momentum of container
                    self._dp_container.append(
                        [np.linalg.norm(self._ball[ball_1]._dp), self._global_time]
                    )
            else:  # Collision with balls
                self._ball[ball_1].collide(self._ball[ball_2])

            if brownian:
                if ball_1 == 0:
                    record = True

        self._N_collisions += 1

        if brownian:
            if record:
                self.record_brownian()
                
    #%% Simulation Randomisation
    ### SIMULATION RANDOMISATION METHODS
    # Gives all the randomisation methods for the simulation.
    def generator_random_position(self, start=0):
        """
        Generates random positions for balls such that they do not overlap.
        
        Parameters:
            start (boolean, optional): The starting index of ball to set random 
                positions for. Used when initialising brownian motion 
                investigation because the large ball starts at the origin.

        Raises:
            Exception: When the balls cannot fit in the container. Reduce
                number of balls or increase container radius.
        """
        for i in range(start, self._N_balls):
            pos = np.zeros(2)
            false_count = 0
            while True:
                if false_count > 1e6:
                    raise Exception("Area of container is too small for ball size")
                x = gen_random_uniform(self._r_container - self._ball[i]._radius)
                y = gen_random_uniform(self._r_container - self._ball[i]._radius)
                while (
                    np.sqrt(x ** 2 + y ** 2)
                    >= self._r_container - self._ball[i]._radius
                ):
                    x = gen_random_uniform(self._r_container - self._ball[i]._radius)
                    y = gen_random_uniform(self._r_container - self._ball[i]._radius)
                pos = np.array([x, y])
                append = False
                for j in range(0, i):
                    distance = np.sqrt(
                        (self._ball[j]._pos[0] - pos[0]) ** 2
                        + (self._ball[j]._pos[1] - pos[1]) ** 2
                    )
                    if distance <= self._ball[i]._radius + self._ball[j]._radius:
                        append = False
                        false_count += 1
                        break
                    else:
                        append = True
                if append or i == 0:
                    break
            self._ball[i].set_pos(pos)

    def generator_random_vel(self, max_speed):
        """
        Generates and sets random velocities for all the balls from a uniform 
            random distribution of x- and y- velocity components.
        
        Parameters:
            max_speed (float): The range of x- and y- velocities component to
                be generated from.
        """
        l = generate_random_vel(self._N_balls, self._random_speed)
        self.set_vel_ball(l)
        
    
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
        self._ball[0].set_pos(np.array([0.0,0.0]))
        self._ball[0].set_radius(radius)
        self._ball[0].set_mass(mass)
        self.generator_random_position(start = 1)


# CHECK
    def set_vel_ball(self, l_vel):
        """
        Sets the velocities of all balls with a given list of velocities.

        Parameters:
            l_vel (list of numpy.ndarray of float): List of the ball velocities
                in their x- and y- directions.
        """
        for i, vel in enumerate(l_vel):
            self._ball[i].set_vel(vel)
# CHECK
    def init_patch(self):
        """
        Initialising the balls and the container patches in the animation.
        Balls and container are drawn using matplotlib.pyplot.Circle objects.
        """
        b_patch = []  # List containing ball patches
        pos_c = self._container._pos
        r_c = self._r_container
        c_outline = plt.Circle(pos_c, r_c, ec="b", fill=False, ls="solid")

        for i, ball in enumerate(self._ball):
            pos_b = ball._pos
            r_b = ball._radius

            if i != 0:  # Generating random colours for patches
                b_patch.append(
                    plt.Circle(
                        pos_b,
                        r_b,
                        ec="black",
                        fc=tuple(
                            (np.random.rand(), np.random.rand(), np.random.rand())
                        ),
                    )
                )

            # Setting first ball to be yellow for visibility in tracing
            # Brownian Motion
            else:
                b_patch.append(plt.Circle(pos_b, r_b, ec="black", fc="yellow"))
        self._b_patch = b_patch
        self._c_outline = c_outline
# CHECK
    def draw(self):
        """
        Drawing the current state of the simulation. Does not animate.
        """
        self.init_patch()

        plt.figure(num="Simulation State")
        ax = plt.axes(
            xlim=(-self._r_container, self._r_container),
            ylim=(-self._r_container, self._r_container),
            aspect="equal",
        )

        ax.add_patch(self._c_outline)  # Drawing container
        for patch in self._b_patch:
            ax.add_patch(patch)  # Drawing balls

        plt.show()
# CHECK

    #%% Simulation Recording ✔️
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
        for i, ball in np.ndenumerate(self._ball):
            KE[i] = 0.5 * ball._mass * np.linalg.norm(ball._vel)**2
        temperature = (np.sum(KE))/(self._N_balls*spc.Boltzmann)
        self._temperature.append([temperature, self._global_time, self._N_collisions])
        
        if len(self._KE) == self._collisions + 1:
            self._KE = pd.DataFrame(self._KE, columns=["T","t","collision"])

    def record_temperature_moyen(self):
        """ record_temperature_moyen | Records the average temperature of the system.
            RECORDS
                (float): Average systemic temperature.
        """
        self._temperature_moyen = np.mean(self._KE["T"]) # Change potentially

    def record_speed(self):
        """ record_speed | Record speed for all balls, for all collisions.
            RECORDS
                (list(float)): Speed of balls, for all collisions.
        """
        for ball in self._ball:
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
        total_KE = np.sum([0.5 * self._m_balls * np.dot(ball._vel,ball._vel) for ball in self._ball])
        self._KE.append([total_KE, self._global_time, self._N_collisions])
        
        if len(self._KE) == self._collisions + 1:
            self._KE = pd.DataFrame(self._KE, columns=['KE','t','collision'])

    def record_distance_absolute(self):
        """
        Writes distances from origin of all balls for all collisions.

        Data Recorded:
            (list of float): Ball distances from the centre of container for
                all collisions.
        """
        for ball in self._ball:
            self._distance_absolute.append(bl.mag_vector(ball._pos))

    def record_distance_relative(self):
        """
        Writes the relative distances between all balls for all collisions.

        Data Recorded:
            (list of float): Relative distances between all pairs of balls for
                all collisions.
        """
        for _, pair in enumerate(self._pairs):
            ball_A = pair[0]
            ball_B = pair[1]
            rel_dist = bl.mag_vector(self._ball[ball_A]._pos - self._ball[ball_B]._pos)
            self._distance_relative.append(rel_dist)

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
        for i, ball in np.ndenumerate (self._ball):
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
                           test_temperature = False,
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
        if temperature or test_temperature:         self.record_temperature()
        if dataset:             self.record_dataset()

    def record_data_pressures(self, pressure=False, test_pressure=False):
        """ record_data_pressures | Records systemic pressure datasets.
            PARAMETERS
                -> pressure(boolean, optional): If True, writes dataset for the average systemic pressure.
                -> test_pressure(boolean, optional): If True, writes dataset of pressure for every _(*)_ container collisions.
        """ 
        if pressure:        self.record_pressure_moyen()
        if test_pressure:   self.record_pressure()

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
        self._pq = hd.heapdict()
        self._collisions = collisions

        print("Starting {self._N_balls} balls, r_balls = {self._r_balls}, \
            speed range = {self._speed_range}, r_container = {self._r_container}, \
            {self._collisions} collisions.")

        if dataset:
            self._dataset = np.zeros((self._N_balls * (self._collisions + 1), 9))

        # Initialising animation
        if animate:
            self.init_patch()

            plt.figure(num="Simulation Animation")
            plt.rcParams.update(plt.rcParamsDefault)
            ax = plt.axes(
                xlim=(-self._r_container, self._r_container),
                ylim=(-self._r_container, self._r_container),
                aspect="equal",
            )

            ax.add_patch(self._c_outline)
            for patch in self._b_patch:
                ax.add_patch(patch)
            plt.pause(time)

        if brownian:
            self.record_brownian()

        self.record_data_states(
            distance_absolute=distance_absolute,
            speed=speed,
            KE=KE,
            test_temperature=test_temperature,
            temperature=temperature,
            distance_relative=distance_relative,
            dataset=dataset,
        )

        # Running first collision
        self.init_collision_time()
        self.init_next_event()
        self.move_balls()

        self._global_time = self._min_dt

        if animate:
            self.update_patch()
            if brownian:
                path = self.trace_brownian()
                ax.add_line(path)
            plt.pause(time)

        self.collide_balls(pressure, test_pressure, brownian)

        self.record_data_states(
            distance_absolute=distance_absolute,
            speed=speed,
            KE=KE,
            test_temperature=test_temperature,
            temperature=temperature,
            distance_relative=distance_relative,
            dataset=dataset,
        )

        if progress_bar:
            self._time_epoch = tm.time()

        for i in range(2, collisions + 1):
            if progress_bar:
                progress(self._time_epoch, i, collisions)
            self.collision_time()
            self.next_event()
            self.move_balls()

            self._global_time = self._min_dt

            if animate:
                self.update_patch()
                if brownian:
                    path = self.trace_brownian()
                    ax.add_line(path)
                plt.pause(time)

            self.collide_balls(pressure, test_pressure, brownian)

            self.record_data_states(
                distance_absolute=distance_absolute,
                speed=speed,
                KE=KE,
                test_temperature=test_temperature,
                temperature=temperature,
                distance_relative=distance_relative,
                dataset=dataset,
            )

        if animate:
            plt.show()

        if temperature:
            self.record_temperature_moyen()

        self.record_data_pressures(pressure=pressure, test_pressure=test_pressure)

        if brownian:
            self.record_brownian(df=True)

        d_output = self.gloss(
            distance_absolute=distance_absolute,
            distance_relative=distance_relative,
            test_pressure=test_pressure,
            speed=speed,
            KE=KE,
            test_temperature=test_temperature,
            temperature=temperature,
            pressure=pressure,
            dataset=dataset,
            brownian=brownian,
        )

        print(
            f"end of {self._N_balls} balls, r_balls = {self._r_balls}, speed range = {self._random_speed}, r_container = {self._r_container}, {self._collisions} collisions"
        )

        return d_output

    #%% Simulation Miscellaneous
    ### SIMULATION MISCELLANEOUS METHODS

def gen_random_uniform(maximum_range):
    """ gen_random_uniform | Provides a uniform distribution centered at (0,0), generating random floats.
        PARAMETERS
            maximum_range (float): sets the maximum range for the distribution.
        RETURNS
            (float): a random float between [-max_range, max_range].
    """
    return np.random.uniform(-maximum_range,maximum_range)

def generate_random_vel(N_balls, random_speed_range):
    """
    Generates random velocities for a given number of balls from a uniform 
    distribution of [-random_speed_range, random_speed_range) for both x- and 
    y- components.

    Parameters:
        N_balls (int): The number of balls.
        random_speed_range (float): The range of speed in the x- and y- 
            velocity component.
    
    Returns:
        l (list of numpy.ndarray of float): List of velocities for all balls.
    """
    l = []
    for _ in range(N_balls):
        l.append(
            np.array([gen_random_uniform(random_speed_range), gen_random_uniform(random_speed_range)])
        )
    return l

def progress(start_time, iterations, all_iterations, desc="Collisions"):
    """ progress | A progress bar to show the progression of an operation.
    PARAMETERS
        start_time
        it
        max_it
        desc
    """    
    current_time = tm.time()

    t_elapsed = current_time - start_time
    t_remaining = t_elapsed * all_iterations / iterations - t_elapsed
    t_elapsed_str = tm.strftime("%H:%M:%S", tm.gmtime(t_elapsed))
    t_remaining_str = tm.strftime("%H:%M:%S", tm.gmtime(t_remaining))

    percentage = round(iterations / all_iterations * 100)
    num_blocks = int(np.floor(percentage / 5))
    blocks = num_blocks * "\u2588" + (20 - num_blocks) * " " # Alter
    desc_str = f"{desc}:"

    # Calculating iterations per second
    if t_elapsed == 0:  it_s = "0.00"
    else:               it_s = round(iterations / t_elapsed, 2)

    progress = f"\r|{blocks}| {percentage}% | {desc_str} {iterations}/{all_iterations} | {t_elapsed_str}/{t_remaining_str} | {it_s}it/s"

    if iterations == all_iterations:    progress += "\n"
    sys.stdout.write(progress)

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


