""" BALL MODULE (import ball as bl) | Created by 𝑀𝒶𝓇𝓉𝒾𝓃 𝒜. 𝐻𝑒, 2022.11.20
    For the simulation of elastic collisions of balls with others, and container.
"""
import numpy as np
from copy import deepcopy

class Ball:
    time = 0
    def __init__(self, mass = 1, radius = 1, pos = np.array([0.0,0.0]),
                 vel = np.array([0.0,0.0]), n_collision = 0, count = 0):
        self._mass = mass                   # Mass of the ball.
        self._radius = radius               # Radius of the ball.
        self._pos = pos                     # Position of the ball, as np.ndarray.
        self._vel = vel                     # Velocity of the ball, as np.ndarray.
        self._n_collision = n_collision     # Number of collisions.
        self._dp = 0                        # Change in momentum of the ball.
        self._count = count                 # Count of the ball.
        return
    
    def __repr__(self):
        return ("BALL PROPERTIES: mass = {self._mass}, radius = {self._radius},\
                 position = {self._pos}, velocity = {self._vel}")

    def __str__(self):
        return ("BALL: m = {self._mass}; r = {self._radius}; x = {self._pos};\
                 v = {self._vel}")
        
    def time_to_collision(self,other):
        r, v = self._pos - other._pos, self._vel - other._vel
        a, b  = np.dot(v,v), 2 * np.dot(r,v)

        # Both balls are stationary (w.r.t. each other).
        if a == 0:  return np.inf
        else:
            if isinstance(other,Container): R = self._radius - other._radius
            else:                           R = self._radius + other._radius

            c = np.dot(r,r) - R**2
            dt = -b/(2*a); discrim = b**2 - 4*a*c
            if not isinstance(other,Container):
                if np.abs(c) <= 10e-13: return np.inf
                    # Encounters a floating point error - no collision.

            if discrim < 0:
                # Complex dt - no collision.
                return np.inf
            elif discrim == 0:              # Repeated root
                if dt <= 0: return np.inf   # Negative, or no dt - no collision.
                else: return dt             # Positive dt.
            elif discrim > 0:               # Two real roots - dt_1, dt_2.
                dt_1 = (-b + np.sqrt(discrim))/(2*a)
                dt_2 = (-b - np.sqrt(discrim))/(2*a)

                if isinstance(other,Container): return np.amax(np.array([dt_1,dt_2]))
                else:
                    if dt_1 <= 0 and dt_2 <= 0: return np.inf
                    elif dt_1 > 0 and dt_2 > 0: return np.amin(np.array([dt_1,dt_2]))
                    elif dt_1 < 0:              return dt_2
                    elif dt_2 < 0:              return dt_1

    def collide(self, other):
        # Transform the positions relative to each other.
        r = self._pos - other._pos
        v = self._vel - other._vel

        # Resolve the velocity into // and |_ to the line of centres.
        u_self_par = np.vdot(self._vel, r) / np.dot(r,r) * r
        u_self_per = self._vel - u_self_par

        # Consider collisions with the container.
        if isinstance(other,Container):
            v_self_par,v_other_par = -u_self_par, np.zeros(2)
            v_self_per = u_self_per
            self.set_vel(v_self_par + v_self_per)
            other.set_vel(np.zeros(2))
            self._count += 1; other._count += 1

            vel_f = self._vel
            vel_i = self._vel
            dv = vel_f - vel_i
            self._dp = dv * self._mass

        # Consider collisions with another ball.
        else:
            m1,m2 = self._mass, other._mass
            # self.set_vel(self._vel - ((2*m2)/(m1*m2))*(np.dot(v,r)/np.dot(r,r))*r)
            # other.set_vel(other._vel - ((2*m1)/(m1*m2))*(np.dot(-v,-r)/np.dot(-r,-r))*(-r))

            # Only care about parallel components in 1D. Consider relative velocities.
            u_other_par = np.vdot(other._vel, r) / np.dot(r,r) * r
            u_other_per = other._vel - u_other_par
            u_self_par_relative = u_self_par - u_other_par
            v_self_par_relative = (m1 - m2)/(m1 + m2) * u_self_par_relative

            v_self_par = v_self_par_relative + u_other_par
            v_other_par = u_self_par - u_other_par + v_self_par
            v_self_per, v_other_per = u_self_per, u_other_per
            self.set_vel(v_self_par + v_self_per)
            other.set_vel(v_other_par + v_other_per)
            self._count += 1; other._count += 1

    def move(self,dt):
        """ move | Move the ball to a new position: r' = r + v * dt.
        """
        self.set_pos(self._pos + dt * self._vel)

    def copy(self):     return deepcopy(self)
    def mass(self):     return self._mass
    def radius(self):   return self._radius
    def pos(self):      return self._pos
    def vel(self):      return self._vel

    def set_mass(self, mass):       self._mass = mass
    def set_radius(self,radius):    self._radius = radius
    def set_pos(self, pos):         self._pos = np.array(pos)
    def set_vel(self,vel):          self._vel = np.array(vel)

class Container(Ball):
    def __init__(self, radius=10, mass=100):
        super().__init__()
        self._radius = radius
        self._mass = mass
        self._count = -1