"""
SIMULATION MODULE (import simulation as sm) | Created by 𝑀𝒶𝓇𝓉𝒾𝓃 𝒜. 𝐻𝑒, 2022.11.20
For the simulation of elastic collisions of balls with others, and container.
"""
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
        self._collision_count = 0
                
        for _ in range(0, N_balls):
            self._ball.append(bl.Ball(radius=r_balls, mass=m_balls))
        if random_position: self.generator_random_position()
        self.generator_random_vel(max_speed=random_speed_range)
        
    ### SIMULATION INFORMATION METHODS
    # Gives all the simulation methods for returning information on the Simulation.
    def __repr__(self):
        return ("Simulation properties: N = {self._N_balls} balls, r_balls = {self._r_balls}, m_balls = {self._m_balls}, r_container = {self._r_container}")
    
    def __str__(self):
        return ("Simulation: N = {self._N_balls} balls, r_balls = {self._r_balls}, m_balls = {self._m_balls}, r_container = {self._r_container}")

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
                x = rng_uniform(self._r_container - self._ball[i]._radius)
                y = rng_uniform(self._r_container - self._ball[i]._radius)
                while (
                    np.sqrt(x ** 2 + y ** 2)
                    >= self._r_container - self._ball[i]._radius
                ):
                    x = rng_uniform(self._r_container - self._ball[i]._radius)
                    y = rng_uniform(self._r_container - self._ball[i]._radius)
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
        self._balls[0].set_pos(np.array([0.0,0.0]))
        self._balls[0].set_radius(radius)
        self._balls[0].set_mass(mass)
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

# CHECK
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

        self._collision_count += 1

        if brownian:
            if record:
                self.record_brownian()
                

    ### SIMULATION RECORDING METHODS
    # Gives simulation methods for recording data.
    def record_dataset(self):
        """
        Writes all simulation information into a pandas.DataFrame.

        Data recorded:
            (pandas.DataFrame of [ball, mass, x, y, vx, vy, collision, t, 
            container]):
                ball (int): Ball number.
                mass (float): Mass of ball.
                x (float): x-coordinate of ball.
                y (float): y-coordinate of ball.
                vx (float): x-velocity of ball.
                vy (float): y-velocity of ball.
                collision (int): Collision number.
                t (float): Time of collision.
                container (boolean): If True, the ball collided with the 
                    container.
        """

        for i, ball in enumerate(self._ball):
            j = self._collision_count * self._N_balls
            k = j + i
            self._dataset[k, 0] = i
            self._dataset[k, 1] = ball._mass
            self._dataset[k, 2] = ball._pos[0]
            self._dataset[k, 3] = ball._pos[1]
            self._dataset[k, 4] = ball._vel[0]
            self._dataset[k, 5] = ball._vel[1]
            self._dataset[k, 6] = self._collision_count
            self._dataset[k, 7] = self._global_time

            # Checks if it is a collision with the container
            if (
                np.abs(bl.mag_vector(ball._pos) - self._r_container + ball._radius)
                <= 10e-10
            ):
                self._dataset[k, 8] = True
            else:
                self._dataset[k, 8] = False

        # Writes the complete data into a pandas.DataFrame
        if self._collision_count == self._collisions:
            self._dataset = pd.DataFrame(
                self._dataset,
                columns=[
                    "ball",
                    "mass",
                    "x",
                    "y",
                    "vx",
                    "vy",
                    "collision",
                    "t",
                    "container",
                ],
            )

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

    def record_speed(self):
        """
        Writes the speeds of all balls for all collisions.

        Data Recorded:
            (list of float): Ball speeds for all collisions.
        """
        for ball in self._ball:
            self._speed.append(bl.mag_vector(ball._vel))

    def record_pressure(self):
        """
        Writes pressure values for every 50 container collisions.

        Raises:
            IndexError: If the number of collisions are insufficient to 
                calculate pressure.

        Data Recorded:
            (pandas.DataFrame of [pressure,t]):
                pressure (float): Pressure of the system.
                t (float): Time.
        """

        self._pressure = []

        N_coll = 50  # Number of collisions to average over

        if not isinstance(self._dp_container, np.ndarray):
            self._dp_container = np.array(self._dp_container)

        try:
            max_t = self._dp_container[-1, 1]
        except IndexError:  # No collisions with the container
            print("Number of collisions insufficient for pressure")
            self._pressure = np.nan
            return

        # Only picks the last 80% of data when system achieves steady state
        min_t = max_t / 5
        start = 0

        # Determining starting index of pressure data
        while self._dp_container[start, 1] <= min_t:
            start += 1
            if start == len(self._dp_container) - 1:
                print("Number of collisions insufficient for pressure")
                self._pressure = np.nan
                return
        start += (len(self._dp_container) - start) % N_coll
        new_dp = self._dp_container[start:, :]

        # Calculating pressure over different times
        N_pressure = int(len(new_dp) / N_coll)
        for i in range(N_pressure):
            index = i * N_coll
            # Takes sum of momentum change and divide by time period
            pressure = np.sum(new_dp[index : index + N_coll - 1, 0]) / (
                (new_dp[index + N_coll - 1, 1] - new_dp[index, 1])
                * 2
                * np.pi
                * self._r_container
            )
            time = (new_dp[index + N_coll - 1, 1] + new_dp[index, 1]) / 2
            self._pressure.append([pressure, time])

        self._pressure = pd.DataFrame(self._pressure, columns=["pressure", "t"])

    def record_pressure_moyen(self):
        """
        Records average pressure of the system.
        Only gives meaningful values if Simulation.run(pressure_moyen=True).

        Raises:
            IndexError: If the number of collisions are insufficient to 
                calculate pressure.

        Data Recorded:
            (float): Average steady state pressure of the system. 
        """
        if not isinstance(self._dp_container, np.ndarray):
            self._dp_container = np.array(self._dp_container)

        try:
            max_t = self._dp_container[-1, 1]
        except IndexError:  # No collision with the container
            print("Number of collisions insufficient for pressure")
            self._pressure_moyen = np.nan
            return

        # Only picks the last 80% of data when system achieves steady state
        min_t = max_t / 5
        start = 0
        # Determining starting index of pressure data
        while self._dp_container[start, 1] <= min_t:
            start += 1
            if start == len(self._dp_container) - 1:
                print("Number of collisions insufficient for pressure")
                self._pressure_moyen = np.nan
                return
        min_t = self._dp_container[start, 1]

        # Average pressure is sum of momentum change divided by time period
        self._pressure_moyen = np.sum(self._dp_container[start:, 0]) / (
            (max_t - min_t) * 2 * np.pi * self._r_container
        )

    def record_KE(self):
        """
        Record total kinetic energy of the system for all collisions.

        Data Recorded:
            (pandas.DataFrame of [KE, t, collision]):
                KE (float): Kinetic Energy of the system.
                t (float): Time.
                collision (int): Collision number.
        """
        KE = np.sum(
            [0.5 * self._m_balls * bl.magsquare_vector(ball._vel) for ball in self._ball]
        )
        self._KE.append([KE, self._global_time, self._collision_count])

        # Writing completed data into pandas.DataFrame
        if len(self._KE) == self._collisions + 1:
            self._KE = pd.DataFrame(self._KE, columns=["KE", "t", "collision"])

    def record_temperature(self):
        """
        Gives temperature of the system at all collision times.

        Data Recorded:
            (pandas.DataFrame of [T, t, collision]):
                T (float): Temperature of the system.
                t (float): Time.
                collision (int): Collision number.
        """
        kb = 1.38064852e-23
        KE = np.zeros(self._N_balls)
        for i, ball in enumerate(self._ball):
            KE[i] = 0.5 * ball._mass * np.linalg.norm(ball._vel) ** 2
        temperature = (np.sum(KE)) / (self._N_balls * kb)

        self._temperature.append(
            [temperature, self._global_time, self._collision_count]
        )

        # Writing the completed data into pandas.DataFrame
        if len(self._temperature) == self._collisions + 1:
            self._temperature = pd.DataFrame(
                self._temperature, columns=["T", "t", "collision"]
            )

    def record_temperature_moyen(self):
        """
        Records average temperature of the system.

        Parameters:
            df (boolean, optional): If True, converts the data into pandas.
                DataFrame. Triggered at the end of collisions

        Data Recorded:
            (float): Average temperature of the system.
        """
        self._temperature_moyen = np.mean(self._temperature["T"])

    def record_brownian(self, df=False):
        """
        Writes dataset required for Brownian Motion investigation.

        Data Recorded:
            (pandas.DataFrame of [x, y, t, collision, hit]):
                x (float): x-coordinate of ball.
                y (float): y-coordinate of ball.
                t (float): Time of collision.
                collision (float): Collision number.
                hit (boolean): True if the collision that took place hit the 
                    ball under investigation.
        """
        if not df:
            self._brownian.append(
                np.array(
                    [
                        self._ball[0]._pos[0],
                        self._ball[0]._pos[1],
                        self._global_time,
                        self._collision_count,
                    ]
                )
            )
        else:  # Writes completed data into pandas.DataFrame
            self._brownian = pd.DataFrame(
                self._brownian, columns=["x", "y", "t", "collision"]
            )

    def record_data_states(
        self,
        distance_absolute=False,
        speed=False,
        KE=False,
        test_temperature=False,
        temperature=False,
        distance_relative=False,
        dataset=False,
    ):
        """
        Writes datasets for the simulation.

        Parameters:
            distance_absolute (boolean, optional): If True, writes dataset for 
                distances of all balls from the origin for all collsions.
            speed (boolean, optional): If True, writes dataset for speeds of 
                all balls for all collisions.
            KE (boolean, optional): If True, writes dataset for kinetic energy 
                of the system for all collision times.
            test_temperature (boolean, optional): If True, writes dataset for 
                temperature of the system at all collision times.
            temperature (boolean, optional): If True, writes dataset for 
                average temperature of the system.
            distance_relative (boolean, optional): If True, writes dataset for relative 
                distances between all balls for all collisions.
            dataset (boolean, optional): If True, writes dataset for all 
                information of the simulation.
        """
        if distance_absolute:
            self.record_distance_absolute()
        if speed:
            self.record_speed()
        if KE:
            self.record_KE()
        if test_temperature or temperature:
            self.record_temperature()
        if distance_relative:
            self.record_distance_relative()
        if dataset:
            self.record_dataset()

    def record_data_pressures(self, pressure=False, test_pressure=False):
        """
        Writes pressure datasets for the system.

        Parameters:
            pressure (boolean, optional): If True, writes dataset for average 
                pressure of the system.
            test_pressure (boolean, optional): If True, writes dataset for 
                pressure of every 100 collisions with the container.
        """
        if pressure:
            self.record_pressure_moyen()

        if test_pressure:
            self.record_pressure()

    def append_data(
        self,
        distance_absolute=False,
        distance_relative=False,
        test_pressure=False,
        speed=False,
        KE=False,
        test_temperature=False,
        temperature=False,
        pressure=False,
        dataset=False,
        brownian=False,
    ):
        """
        Appends required data into a dictionary to be returned at the end of 
        the simulation.

        Parameters:
            distance_absolute (boolean, optional): If True, appends to dictionary 
                distance to centre dataset.
            speed (boolean, optional): If True, appends to dictionary speeds 
                dataset.
            KE (boolean, optional): If True, appends to dictionary kinetic 
                energy dataset.
            test_temperature (boolean, optional): If True, appends to 
                dictionary temperature dataset.
            temperature (boolean, optional): If True, appends to dictionary 
                average temeprature dataset.
            distance_relative (boolean, optional): If True, appends to dictionary 
                relative distances dataset.
            dataset (boolean, optional): If True, appends to dictionary dataset 
                for all information of the simulation.
            brownian (boolean, optional): If True, appends to dictionary 
                dataset for Brownian Motion investigation.

        """
        d_output = {}

        if distance_absolute:
            d_output["distance from centre"] = self._distance_absolute
        if distance_relative:
            d_output["relative distance"] = self._distance_relative
        if test_pressure:
            d_output["pressure"] = self._pressure
        if speed:
            d_output["speed"] = self._speed
        if KE:
            d_output["KE"] = self._KE
        if test_temperature:
            d_output["temperature"] = self._temperature
        if temperature:
            d_output["average temperature"] = self._temperature_moyen
        if pressure:
            d_output["average pressure"] = self._pressure_moyen
        if dataset:
            d_output["dataset"] = self._dataset
        if brownian:
            d_output["brownian"] = self._brownian

        return d_output

    titleFont =     {'fontname': 'Kinnari', 'size': 13}
    axesFont =      {'fontname': 'Kinnari', 'size': 9}
    ticksFont =     {'fontname': 'SF Mono', 'size': 7}
    errorStyle =    {'mew': 1, 'ms': 3, 'capsize': 3, 'color': 'blue', 'ls': ''}
    pointStyle =    {'mew': 1, 'ms': 3, 'color': 'blue'}
    lineStyle =     {'linewidth': 0.5}
    lineStyleBold = {'linewidth': 1}
    histStyle =     {'facecolor': 'green', 'alpha': 0.5, 'edgecolor': 'black'}
    
    ### SIMULATION RUN METHOD
    # The method to run simulations. 
    def run(
        self,
        collisions=10,
        time=0.001,
        *,
        animate=False,
        distance_absolute=False,
        distance_relative=False,
        speed=False,
        pressure=False,
        dataset=False,
        test_pressure=False,
        KE=False,
        test_temperature=False,
        temperature=False,
        brownian=False,
        progress=True,
    ):
        """
        Runs the 2D rigid disc particle collision simulation.


        Parameters:
            collisions (int, optional): Number of collisions in the simulation.
            time (float, optional): Time period between animation frames.
            animate (boolean, optional): If True, produces animation. 
            distance_absolute (boolean, optional): If True, records distances of 
                all balls from origin for all collisions.
            distance_relative (boolean, optional): If True, records relative distances 
                between all balls for all collisions.
            speed (boolean, optional): If True, records speeds of all balls for 
                all collisions.
            pressure (boolean, optional): If True, records pressure for every 
                100 collisions with container.
            dataset (boolean, optional): If True, records dataset of all 
                information of the simulation.
            test_pressure (boolean, optional): If True, records average 
                pressure of the system.
            KE (boolean, optional): If True, records kinetic energy of the 
                system for every collision.
            test_temperature (boolean, optional): If True, records temperature 
                of the system for every collision.
            temperature (boolean, optional): If True, records average 
                temperature of the system.
            brownian (boolean, optional): If True, records data for Brownian 
                Motion investigation.
            progress (boolean, optional): If True, displays the progress bar.

        Returns:
            d_output (dict): Dictionary that contains the required datasets.
        """
        self._global_time = 0
        self._pq = hd.heapdict()
        self._collisions = collisions

        print(
            f"starting {self._N_balls} balls, r_balls = {self._r_balls}, speed range = {self._random_speed}, r_container = {self._r_container}, {self._collisions} collisions"
        )

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

        if progress:
            self._time_epoch = tm.time()

        for i in range(2, collisions + 1):
            if progress:
                progress_bar(self._time_epoch, i, collisions)
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

        d_output = self.append_data(
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


def progress_bar(start_time, it, max_it, desc="Collisions"):
    """
    A custom progress bar to show the progress of any iteration.

    Parameters:
        start_time (float): The time since epoch given by time.time().
        it (int): Current iteration number.
        max_it (int): Maximum iteration number.
        desc (string): Description of the iteration.
    """
    current_time = tm.time()
    t_elapsed = current_time - start_time
    t_remaining = t_elapsed * max_it / it - t_elapsed

    # Converting time into HH:MM:SS
    t_elapsed_str = tm.strftime("%H:%M:%S", tm.gmtime(t_elapsed))
    t_remaining_str = tm.strftime("%H:%M:%S", tm.gmtime(t_remaining))

    percentage = round(it / max_it * 100)

    # Number of blocks in progress bar to be displayed.
    num_blocks = int(np.floor(percentage / 10))
    blocks = num_blocks * "\u2588" + (10 - num_blocks) * " "

    desc_str = f"{desc}:"

    if t_elapsed == 0:
        it_s = "0.00"
    else:  # Calculating iterations per second
        it_s = round(it / t_elapsed, 2)

    progress = f"\r|{blocks}| {percentage}% | {desc_str} {it}/{max_it} | {t_elapsed_str}/{t_remaining_str} | {it_s}it/s"

    if it == max_it:
        progress += "\n"
    sys.stdout.write(progress)


def rng_uniform(max_range):
    """
    Generates random float values given from a uniform distrubition centred at 
        0.
    
    Parameters:
        max_range (float): The maximum range of the uniform distribution.
    
    Returns:
        (float): A random float from [-max_range, max_range]

    """
    return random.uniform(-max_range, max_range)
# CHECK

def temperature_from_rms_speed(rms_speed, m):
    """
    Calculates temperature from the root-mean-squared speed of particles.

    Parameters:
        rms_speed (float): Root-mean-squared speed of the particles.
        m (float): Mass of the particles.
    
    Returns:
        (float): The temperature of the system.
    """
    k = 1.38064852e-23
    return 0.5 * m * rms_speed ** 2 / k


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
            np.array([rng_uniform(random_speed_range), rng_uniform(random_speed_range)])
        )
    return l


def equal_sampling_brownian(brownian_data, N_samples):
    """    
    Samples the position of a ball at regular time intervals.
    
    Parameters:
        brownian_data (pandas.DataFrame of [x, y, t]): All position data of a 
            ball at times when it experiences a collision.
                x (float): The x-coordinates of the ball.
                y (float): The y-coordinates of the ball.
                t (float): The time of collision.
        N_samples (int): The number of samples to be obtained from the data.
    
    Returns:
        (pandas.DataFrame of [x, y, t]): The position data of a ball at regular 
            time intervals.
                x (float): The x-coordinates of the ball.
                y (float): The y-coordinates of the ball.
                t (float): The time (equally sampled).
    """
    x_pos = np.array(brownian_data["x"])
    y_pos = np.array(brownian_data["y"])
    time = np.array(brownian_data["t"])

    x_samples = np.zeros(N_samples)
    y_samples = np.zeros(N_samples)
    t_samples = np.zeros(N_samples)

    dt = (time[-1] - time[0]) / (N_samples + 1)

    index = 0
    i = 0

    # Sampling the positions using time interval dt
    while i < N_samples:
        t_sampling = time[0] + i * dt

        if t_sampling == time[index]:
            x_samples[i] = x_pos[index]
            y_samples[i] = y_pos[index]
            t_samples[i] = time[index]
            i += 1
        else:
            # Time is larger than next indexed value
            if t_sampling >= time[index + 1]:
                index += 1
            # Time is between the two indexed values
            else:
                t1 = time[index]
                t2 = time[index + 1]
                x1 = x_pos[index]
                x2 = x_pos[index + 1]
                y1 = y_pos[index]
                y2 = y_pos[index + 1]

                # Using a linear equation to interpolate
                x = linear_interpolation(t_sampling, t1, t2, x1, x2)
                y = linear_interpolation(t_sampling, t1, t2, y1, y2)
                x_samples[i] = x
                y_samples[i] = y
                t_samples[i] = t_sampling
                i += 1

    df = pd.DataFrame([x_samples, y_samples, t_samples])
    df = df.transpose()
    df = df.rename(columns={0: "x", 1: "y", 2: "t"})
    return df


def linear_interpolation(x, x1, x2, y1, y2):
    """
    Returns a y-value given a linear equation.

    Parameters:
        x (float): The x-value for the returned result.
        x1 (float): The x-value of known point 1.
        x2 (float): The x-value of known point 2.
        y1 (float): The y-value of known point 1.
        y2 (float): The y-value of known point 2.

    Returns:
        (float): The y-value at point x in the linear equation.
    """
    return y1 + (y2 - y1) / (x2 - x1) * (x - x1)


def trace_paths_brownian(data):
    """
    Traces out the path of a ball using given position data.

    Parameters:
        data (pandas.DataFrame of [x,y]): The time ordered position data of a 
            ball at times of collision.
                x (float): The x-coordinate of the ball.
                y (float): The y-coordinate of the ball.
    
    Returns:
        (list of matplotlib.pyplot.Line2D): The list of Line2D objects that can 
            be plotted on a graph using matplotlib.pyplot
    """
    x = data["x"]
    y = data["y"]
    l_path = []

    for i in range(1, len(x)):
        path = plt.Line2D(
            xdata=[x[i], x[i - 1]], ydata=[y[i], y[i - 1]], color="0.01", alpha=0.2,
        )
        l_path.append(path)

    return l_path



