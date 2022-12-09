    def randomiser(self, start=0, position=False, velocity=False, max_speed=5,
                   N_balls=50):
        """ randomiser | Generates random velocities and positions.
                < PARAMETERS >
                -> start(float, optional): The starting position.
                   N.B. Must have position=True enabled.
                -> position(boolean): Generate the non-overlapping randomised
                                      ball positions.
                -> velocity(boolean): Generate random velocities from a uniform
                        distribution for a given number of balls in the range of:
                        [-maximum_speed_range, maximum_speed_range].
                -> max_speed(float, optional): Give the maxmimum speed generated.                   
                RAISES
                    Exception: Balls cannot fit in this container. Reduce
                               N_balls or increase r_container.
        """
        if velocity:
            list = []
            for _ in range(N_balls):
                list.append(np.array([gen_random_uniform(-max_speed),
                                    gen_random_uniform(max_speed)]))
            for i, vel in enumerate(list):     self._balls[i].set_vel(vel)
        
        if position:
            for i in range(start, N_balls):
                pos = np.zeros(2)
                false_count = 0
                while True:
                    if false_count > 1e6:
                        raise Exception("Area of container is too small for ball size")
                    x = gen_random_uniform(self._r_container - self._balls[i]._radius)
                    y = gen_random_uniform(self._r_container - self._balls[i]._radius)
                    while (
                        np.sqrt(x ** 2 + y ** 2)
                        >= self._r_container - self._balls[i]._radius
                    ):
                        x = gen_random_uniform(self._r_container - self._balls[i]._radius)
                        y = gen_random_uniform(self._r_container - self._balls[i]._radius)
                    pos = np.array([x, y])
                    append = False
                    for j in range(0, i):
                        distance = np.sqrt((self._balls[j]._pos[0] - pos[0])**2
                            + (self._balls[j]._pos[1] - pos[1]) ** 2)
                        if distance <= self._balls[i]._radius + self._balls[j]._radius:
                            append = False
                            false_count += 1
                            break
                        else:
                            append = True
                    if append or i == 0:
                        break
                self._balls[i].set_pos(pos)