"""
BALL MODULE (import ball as bl)
Created by Martin A. He, 2022.11.20
For the simulation of elastic collisions of balls with others, and container.
"""

import numpy as np
from copy import deepcopy

class Ball:
    time = 0
    def __init__(
        self,
        mass = 1,
        radius = 1,
        pos_ball = np.array([0.0,0.0]),
        vel_ball = np.array([0.0,0.0]),
        n_collision = 0,
        count = 0,
    ):
        self._mass = mass
        self._radius = radius
        self._pos_ball = pos_ball
        self._vel_ball = vel_ball
        self._n_collision = n_collision
        self._dp = 0
        self._count = count
        return
    
    ### BALL INFORMATION METHODS
    
    def __repr__(self):
        return ("Ball properties: mass = {self._mass}, radius = {self._radius}, position = {self._pos_ball}, velocity = {self._vel_ball}")
    def __str__(self):
        return ("Ball: m = {self._mass}; r = {self._radius}; x = {self._pos_ball}; v = {self._vel_ball}")

    ### BALL MOVEMENT METHODS
    
    def time_to_collision(self,other):
        """
        Calculate the time until the next collision between this ball and another one, or the container.
        RETURNS
            time (float): the time to the next collision for a ball.
            
        VARIABLES
            R: distance between centre of ball, and centre of other object during collision.
            a, b, c: terms in the quadratic equation for calculation of dt.
            discrim: calculation of the discriminant for dt.
            
        CLASSIFICATION OF SOLUTIONS TO dt
        -> No solution exists for Δ < 0.
        -> A unique solution exists for Δ = 0.
            - If the value of dt > 0, dt is a valid solution.
            - If the value of dt <= 0, dt does not exist.
        -> Two unique solutions exist for Δ > 0.
            - If both dt values are > 0, then the next collision is the smallest of the two.
            - If both dt values are <= 0, dt does not exist.
            - If dt_1 <= 0, but dt_2 > 0, return dt_2.
            - If dt_2 <= 0, but dt_1 > 0, return dt_1.
        """
        r = self._pos_ball - other._pos_ball
        v = self._vel_ball - other._vel_ball
        
        a = np.dot(v,v)
        b = 2 * np.dot(r,v)
        
        if a == 0:  
            # Both balls are stationary (w.r.t. each other).
            return np.nan
        else:
            if isinstance(other,Container):
                R = self._radius - other._radius
            else:
                R = self._radius + other._radius
            c = np.dot(r,r) - R**2
            
            if not isinstance(other,Container):
                if np.abs(c) <= 10e-13:
                    # Encounters a floating point error - no collision.
                    return np.nan
                
            discrim = b**2 - 4*a*c
            if discrim < 0: 
                # Complex dt - no collision.
                return np.nan
            elif discrim == 0:
                # Repeated root
                dt = -b/(2*a)
                if dt <= 0:
                    # Negative dt - no collision.
                    return np.nan
                else:
                    return dt
            elif discrim > 0: 
                # Two real roots - dt_1, dt_2.
                dt_1 = (-b + np.sqrt(discrim))/(2*a)
                dt_2 = (-b - np.sqrt(discrim))/(2*a)
                
                if isinstance(other,Container):
                    return np.amax(np.array([dt_1,dt_2]))
                else:
                    if dt_1 <= 0 and dt_2 <= 0:
                        return np.nan
                    elif dt_1 > 0 and dt_2 > 0:
                        return np.amin(np.array([dt_1,dt_2]))
                    elif dt_1 < 0:
                        return dt_2
                    elif dt_2 < 0:
                        return dt_1

    def collide(self, other):
        """
        Make the changes to the velocities of the ball, and the other one due to a collision.
        PARAMETERS
            other (Ball): the other ball which collides with this one.
        """
        # Transform the positions relative to each other.
        r = self._pos_ball - other._pos_ball
        
        # Resolve the velocity into // and |_ to the line of centres.
        u_self_par = np.vdot(self._vel_ball, r) / np.dot(r,r) * r
        u_self_per = self._vel_ball - u_self_par
        
        # Consider collisions with the container.
        if isinstance(other,Container):
            v_self_par,v_other_par = -u_self_par, np.zeros(2)
            v_self_per = u_self_per
            self.set_vel(v_self_par + v_self_per)
            other.set_vel(np.zeros(2))
            self.count, other.count += 1, 1
            
            vel_f = self._vel_ball
            vel_i = self._vel_ball
            dv = vel_f - vel_i
            self._dp = dv * self._mass
            
        # Consider collisions with another ball.
        else:
            m1,m2 = self._mass, other._mass
            
            # Only care about parallel components in 1D. Consider relative velocities.
            u_other_par,u_other_per = np.vdot(other._vel_ball, r) / np.dot(r,r) * r, other._v_ball - u_other_par
            u_self_par_relative = u_self_par - u_other_par
            v_self_par_relative = (m1 - m2)/(m1 + m2) * u_self_par_relative
            
            v_self_par = v_self_par_relative + u_other_par
            v_other_par = u_self_par - u_other_par + v_self_par
            v_self_per, v_other_per = u_self_per, u_other_per
            self.set_vel(v_self_par + v_self_per)
            other.set_vel(v_other_par + v_other_per)
            self._count, other._count += 1, 1
            
    def move(self,dt):
        """
        Move the ball to a new position: r' = r + v * dt.
        """
        self.set_pos(self._pos + dt * self._vel)


    ### BALL ATTRIBUTE METHODS
    
    def mass(self):
        """
        Return the mass of the ball.
        RETURNS
            mass (float): the mass of the ball.
        """
        return self._mass
    
    def radius(self):
        """
        Return the radius of the ball.
        RETURNS
            radius (float): the radius of the ball.
        """
        return self._radius
    
    def pos(self):
        """
        Return the current position.
        RETURNS 
            position (np.ndarray w/ float values): the current position of the centre of the ball.
        """
        return self._pos_ball

    def vel(self):
        """
        Return the current velocity.
        RETURNS
            velocity (np.ndarray w/ float values): the current 2D velocity of the ball.
        """
        return self._vel_ball
    
    def set_mass(self, mass):
        """
        Define the new mass of the ball.
        PARAMETERS
            mass (float): the new defined mass of the ball.
        """
        self._mass = mass
    
    def set_radius(self,radius):
        """
        Define the new radius of the ball.
        PARAMETERS
            radius (float): the new defined radius of the ball.
        """
        self._radius = radius
        
    def set_pos(self, pos):
        """
        Define the new position of the ball.
        PARAMETERS
            pos (np.ndarray(contains: float)): the new defined position of the ball.
        """
        self._pos_ball = np.array(pos)
    
    def set_vel(self,vel):
        """
        Define the new velocity of the ball.
        PARAMETERS
        
        """
        self._vel_ball = np.array(vel)
        
    ### BALL OTHER METHODS
    
    def copy(self):
        """
        Performs a deepcopy of the ball.
        
        RETURNS
            (Ball): an identical Ball object, located within itself.
        """
        return deepcopy(self)

# Inherit the ball class for the container class.
class Container(Ball):
    """
    This is a Container class. The container is used to enclose the balls in
    the 2D rigid disc collision.

    PARAMETERS
        radius (float): the radius of the container.
        mass (float): the mass of the container.
    """

    def __init__(self, radius=10, mass=100):
        super().__init__()
        self._radius = radius
        self._mass = mass
        self._count = -1