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
                
    def plot_van_der_waals_law(m_ball=5e-26, r_container=5, N_balls=300,
                               collisions = 10000, r_balls=np.linspace(500,2000,7),
                               random_speed_ranges=np.linspace(500,2000,7)):
        """ plot_van_der_waals_law | Plots the Van der Waals graph.
        """
        def van_der_waals_law(T, N, V, b):
            """ van_der_waals_law | Returns the value from the Law.
                    < PARAMETERS >
                    -> T(float): The system temperature.
                    -> V(float): The container volume.
                    -> N(int): The number of gas particles.
                    -> b(float): The effective area of gas particles.
                    RETURNS
                        (float): The pressure on the container walls.
            """
            return (N.spc.Boltzmann * N)/(V-N*b)*T
        
        def calculate_b(volume_container, N_balls, m):
            """ calculate_b | Calculates the value of effective gas area.
                    < PARAMETERS >
                    -> volume_container(float): The container volume.
                    -> N_balls(int): The number of balls in this system.
                    -> m(float): Gradient of P against T.z
                    RETURNS
                        (float): The effective volume of a gas particle, b.
            """
            return(volume_container/N_balls) - spc.Boltzmann/m
        
        def calculate_unc_b(m, unc_m):
            """ calculate_unc_b | Calculate the uncertainty in b.
                    < PARAMETERS >
            """
            return (spc.Boltzmann/m**2) * unc_m
        
        def power_law(x, n, A, B):
            """ power_law | Calculates a power law.
            """
            return A*x**n + B
        
        def run(params):
            """ run | Runs the simulation.
            """
            r_ball, random_speed_range = params[0], params[1]
            simulation_van_der_waals = Simulation(N_balls = N_balls, r_container=r_container,
                                                  r_balls=r_ball,
                                                  random_speed_range=random_speed_range)
            sim = simulation_van_der_waals.run(collisions=collisions, pressure=True,
                                               temperature=True, progress_bar=False)
            return sim
        
        parameters, pressures, temperatures = [], [], []
        volume = r_container**2 * np.pi
        
        if __name__ == "__main__":
            for r_ball in r_balls:
                for random_speed_range in random_speed_ranges:
                    parameters.append([r_ball,random_speed_range])
                    
            t_start = tm.perf_counter()
            with ft.ProcessPoolExecutor() as executor:
                results = executor.map(run, parameters)
            t_end = tm.perf_counter()
            
            print(f"Time elapsed: {round(t_end-t_start,1)}s")
            for result in results:
                pressures.append(result["average pressure"])
                temperatures.append(result["average temperature"])
                
            list_b, list_V = np.zeros(len(r_balls)),  np.zeros(len(r_balls))
            unc_b = np.zeros(len(r_balls))
            
            for i, r_ball in enumerate(r_balls):        # Data to plot b vs. V.
                temperatures_cache, pressures_cache = [], []
                for j, random_speed_range in enumerate(random_speed_ranges):
                    temperatures_cache.append(temperatures[i*len(random_speed_ranges) + j])
                    pressures_cache.append(pressures[i*len(random_speed_ranges) + j])
                linear = sp.stats.linregress(temperatures_cache, pressures_cache)
                
                m, unc_m = linear.slope, linear.stderr  # Linear regression.
                list_b[i], list_V[i] = calculate_b(volume, N_balls, m),\
                                       np.pi * r_ball**2
                unc_b[i] = calculate_unc_b(m, unc_m)
                
            guess_n, guess_A, guess_B = 0.5, 1, 0.5
            p0_power_law = [guess_n, guess_A, guess_B]
            fit_power_law = opt.curve_fit(power_law, list_V, list_b,
                                            p0=p0_power_law, sigma=unc_b)
            
            try:
                print(f"Power of b: {fit_power_law[0][0]} ± {np.sqrt(fit_power_law[1][0,0])}")
                print(f"Scale Factor A: {fit_power_law[0][1]} ± {np.sqrt(fit_power_law[1][1,1])}")
                print(f"Y-Shift B: {fit_power_law[0][2]} ± {np.sqrt(fit_power_law[1][2,2])}")
            except ZeroDivisionError:
                print("Encounter a Zero Division Error.")
            
            array_fit = np.linspace(np.amin(list_V), np.amax(list_V), 1000)
            data_fit_power_law = power_law(array_fit, *fit_power_law[0])
            
            legend_fit = r"Power Law: $n = %s \pm %s$" \
                %(float("%.2g" % fit_power_law[0][0]),
                  float("%.1g" % np.sqrt(fit_power_law[1][0,0])))

            print("Plotting Graph 1 of 1")

            plt.figure(num="Power Law of b and V")
            sns.set(context="paper", style="darkgrid", palette="muted")

            plt.plot(array_fit, data_fit_power_law, label=legend_fit, lw=2, alpha=0.8)
            plt.plot(list_V, list_b, "o", mew=0.5, mec="white")
            plt.errorbar(list_V, list_b, yerr=unc_b, fmt="none", color="black", capsize=3)

            plt.title("Power Law Scaling of Ball Effective Area")
            plt.xlabel(r"Area of 1 Ball $V_{ball}$ /$m^2$")
            plt.ylabel(r"Effective Area of 1 Ball $b$ /$m^2$")
            plt.legend()
            plt.tight_layout()
            plt.show()
