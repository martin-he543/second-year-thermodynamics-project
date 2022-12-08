import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

class Brownian:
    def trace_paths_brownian(data):
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
        def linear_inter(x, x_1, x_2, y_1, y_2):
            """ linear_inter | Returns y-value for two given linear equations.
            """
            return y_1 + (y_2 - y_1) / (x_2 - x_1) * (x - x_1)
        
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
                
                x, y = linear_inter(t_samp,t1,t2,x1,x2), linear_inter(t_samp,t1,t2,y1,y2)
                x_samp[i], y_samp[i], t_samp[i] = x, y, t_samp
                i += 1
                
        return pd.DataFrame([x_samp, y_samp, t_samp]).transpose().rename(columns=\
            {0: "x", 1: "y", 2: "t"})
        
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