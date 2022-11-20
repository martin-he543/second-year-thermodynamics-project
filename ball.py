import numpy as np
from copy import deepcopy

class Ball:
    time = 0
    def __init__(
        self,
        m_ball=1,
        r_ball=1,
        pos_ball=np.array([0.0,0.0]),
        vel_ball=np.array([0.0,0.0]),
        n_collision=0,
        count=0,
    ):
        self._m_ball = m_ball
        self._r_ball = r_ball
        self._pos_ball = pos_ball
        self._vel_ball = vel_ball
        self._n_collision = n_collision
        self._dp = 0
        self._count = count
        return
    
    ### BALL INFORMATION METHODS
    
    def __repr__(self):
        return ("Ball properties: mass = {self._m_ball}, radius = {self._r_ball}, position = {self._pos_ball}, velocity = {self._vel_ball}")
    def __str__(self):
        return ("Ball: m = {self._m_ball}; r = {self._r_ball}; x = {self._pos_ball}; v = {self._vel_ball}")

    ### BALL MOVEMENT METHODS
    
    def time_to_collision(self,other):
        """
        Calculate the time until the next collision between this ball and another one, or the container.
        RETURNS
            time (float): the time to the next collision for a ball.
        """
        r = self._pos_ball - other._pos_ball
        v = self._vel_ball - other._vel_ball
        
        a = np.dot(v,v)
        b = 2 * np.dot(r,v)
        
        if a == 0:  
            # Both balls are stationary (w.r.t. each other).
            return np.inf
        else:
            if isinstance(other,Container):
                R = self._r_ball - other._r_ball
            else:
                R = self._r_ball + other._r_ball
            c = np.dot(r,r) - R**2
            
            if not isinstance(other,Container):
                if np.abs(c) <= 10e-13:
                    # Encounters a floating point error - no collision.
                    return np.inf
                
            discrim = b**2 - 4*a*c
            if discrim < 0: 
                # Complex dt - no collision.
                return np.inf
            elif discrim == 0:
                # Repeated root
                dt = -b/(2*a)
                if dt <= 0:
                    # Negative dt - no collision.
                    return np.inf
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
                        return np.inf
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
        # Transform the velocity to C.O.M. frame.
        r = self._pos_ball - other._pos_ball
        
        # Resolve the velocity into // and |_ to line of centres.
        u1_par = np.vdot(self._vel_ball, r) / np.dot(r,r) * r
        u1_per = self._vel_ball - u1_par
        
        # Consider collisions with the container.
        if isinstance(other,Container):
            v1_par,v2_par = -u1_par, np.zeros(2)
            v1_per = u1_per
            self.set_vel(v1_par + v1_per)
            other.set_vel(np.zeros(2))
            self.count += 1
            other.count += 1
            
            vel_f = self._vel_ball
            vel_i = self._vel_ball
            dv = vel_f - vel_i
            self._dp = dv * self._m_ball
            
        # Consider collisions with another ball.
        else:
            m1,m2 = self._m_ball, other._m_ball
            
            u2_par,u2_per = np.vdot(other._vel_ball, r) / np.dot(r,r) * r, other._v_ball - u2_par
            u1_par_translated = u1_par - u2_par
            v1_par_translated = (m1 - m2)/(m1 + m2) * u1_par_translated
            
            v1_par = v1_par_translated + u2_par
            v2_par = u1_par - u2_par + v1_par
            v1_per = u1_per
            v2_per = u2_per
            self.set_vel(v1_par + v1_per)
            other.set_vel(v2_par + v2_per)
            self._count += 1
            other._count += 1
            
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
            mass (float): The mass of the ball.
        """
        return self._m_ball
    
    def radius(self):
        """
        Return the radius of the ball.
        RETURNS
            radius (float): The radius of the ball.
        """
        return self._r_ball
    
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
        self._m_ball = mass
    
    def set_radius(self,radius):
        """
        Define the new radius of the ball.
        PARAMETERS
            radius (float): the new defined radius of the ball.
        """
        self._r_ball = radius
        
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


class Container(Ball):
    """
    This is a Container class. The container is used to enclose the balls in
    the 2D rigid disc collision.

    PARAMETERS
        radius (float): The radius of the container.
        mass (float): The mass of the container.
    """

    def __init__(self, radius=10, mass=100):
        super().__init__()
        self._radius = radius
        self._mass = mass
        self._count = -1