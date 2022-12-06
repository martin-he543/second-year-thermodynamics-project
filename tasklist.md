Change everything to np array.
Plagiarism code checker.
Switch to np methods.
Create a minified version of the code for running.
Difference between collisions and N_collisions.
Change orders of function lists.


Pressure
Temperature
Speed
KE
Distance
Brownian
Dataset
Data States

        self._balls = []                    # List of container balls.
        self._N_balls = N_balls             # The number of container balls.
        self._N_collisions = 0              # The number of ball collisions.
        self._m_balls = m_balls             # The mass of the container balls.
        self._r_balls = r_balls             # The radius of the container balls.
        self._r_container = r_container     # The radius of the container.
        self._pairs = self.pair_combn()     # List of ball pair combinations.
    # ℭ𝔬𝔩𝔩𝔦𝔰𝔦𝔬𝔫 𝔓𝔯𝔬𝔭𝔢𝔯𝔱𝔦𝔢𝔰:
        self._temperature = []              # System T, ∀ collisions.
        self._KE = []                       # System K.E, ∀ collisions.
        self._speed = []                    # Speed of all balls ∀ collisions.
        self._distance_absolute = []        # Ball distance from origin.
        self._distance_relative = []        # Relative distance between balls.
    # 𝔅𝔞𝔩𝔩 ℜ𝔞𝔫𝔡𝔬𝔪𝔦𝔰𝔞𝔱𝔦𝔬𝔫:
        self._random_position = random_position     # Sets random positioning.
        self._random_speed = random_speed_range     # Sets random speed range.
        self._brownian = []                 # Brownian investigation datasets.
    # ℭ𝔬𝔫𝔱𝔞𝔦𝔫𝔢𝔯 𝔓𝔯𝔬𝔭𝔢𝔯𝔱𝔦𝔢𝔰:
        self._container = bl.Container(radius=r_container) # Chosen container.
        self._dp_container = []             # Changes in container momentum.
        self._N_container_collisions = 0    # Number of container collisions.
    # 𝔗𝔦𝔪𝔢 𝔓𝔯𝔬𝔭𝔢𝔯𝔱𝔦𝔢𝔰:
        self._utc = 0                       # Universal Co-ordinated Time.
        self._min_dt = 0                    # Minimum time to next collision.
    # 𝔈𝔳𝔢𝔫𝔱 𝔓𝔯𝔬𝔭𝔢𝔯𝔱𝔦𝔢𝔰:     
        self._pq = hd.heapdict()            # Priority queue.
        self._events = []                   # List of current events.
  