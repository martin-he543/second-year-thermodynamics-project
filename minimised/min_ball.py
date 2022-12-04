# BALL CLASS, MINIFIED BY ABOUT 70%
import numpy as np; from copy import deepcopy
class Ball:
    time = 0
    def __init__(self,mass = 1,radius = 1,pos_ball = np.array([0.0,0.0]),vel_ball = np.array([0.0,0.0]),n_collision = 0,count = 0):
        self._mass = mass; self._radius = radius; self._pos_ball = pos_ball; self._vel_ball = vel_ball; self._n_collision = n_collision; self._dp = 0; self._count = count
        return
    def __repr__(self): return ("Ball properties: mass = {self._mass}, radius = {self._radius}, position = {self._pos_ball}, velocity = {self._vel_ball}")
    def __str__(self): return ("Ball: m = {self._mass}; r = {self._radius}; x = {self._pos_ball}; v = {self._vel_ball}")    
    def time_to_collision(self,other):
        r,v = self._pos_ball - other._pos_ball,self._vel_ball - other._vel_ball
        a,b = np.dot(v,v),2 * np.dot(r,v); c = np.dot(r,r) - R**2; discrim = b**2 - 4*a*c; dt = -b/(2*a)
        if a == 0: return np.inf
        else:
            if isinstance(other,Container): R = self._radius - other._radius
            else: R = self._radius + other._radius
            if not isinstance(other,Container):
                if np.abs(c) <= 10e-13: return np.inf           
            if discrim < 0: return np.inf
            elif discrim == 0:
                if dt <= 0: return np.inf
                else: return dt
            elif discrim > 0:
                dt_1,dt_2 = (-b + np.sqrt(discrim))/(2*a),(-b - np.sqrt(discrim))/(2*a)
                if isinstance(other,Container): return np.amax(np.array([dt_1,dt_2]))
                else:
                    if dt_1 <= 0 and dt_2 <= 0: return np.inf
                    elif dt_1 > 0 and dt_2 > 0: return np.amin(np.array([dt_1,dt_2]))
                    elif dt_1 < 0: return dt_2
                    elif dt_2 < 0: return dt_1
    def collide(self, other):
        r = self._pos_ball - other._pos_ball; u_self_par = np.vdot(self._vel_ball, r) / np.dot(r,r) * r; u_self_per = self._vel_ball - u_self_par
        if isinstance(other,Container):
            v_self_par,v_other_par = -u_self_par, np.zeros(2); v_self_per = u_self_per; self.set_vel(v_self_par + v_self_per); other.set_vel(np.zeros(2)); self.count, other.count += 1, 1
            vel_f,vel_i = self._vel_ball,self._vel_ball; dv = vel_f - vel_i; self._dp = dv * self._mass
        else:
            m1,m2 = self._mass, other._mass; u_other_par,u_other_per = np.vdot(other._vel_ball, r) / np.dot(r,r) * r, other._v_ball - u_other_par
            u_self_par_relative = u_self_par - u_other_par; v_self_par_relative = (m1 - m2)/(m1 + m2) * u_self_par_relative
            v_self_par = v_self_par_relative + u_other_par; v_other_par = u_self_par - u_other_par + v_self_par; v_self_per, v_other_per = u_self_per, u_other_per
            self.set_vel(v_self_par + v_self_per); other.set_vel(v_other_par + v_other_per); self._count, other._count += 1, 1
    def move(self,dt): self.set_pos(self._pos + dt * self._vel)
    def mass(self): return self._mass
    def radius(self): return self._radius
    def pos(self): return self._pos_ball
    def vel(self): return self._vel_ball
    def set_mass(self, mass): self._mass = mass
    def set_radius(self,radius): self._radius = radius
    def set_pos(self, pos): self._pos_ball = np.array(pos)
    def set_vel(self,vel): self._vel_ball = np.array(vel)
    def copy(self): return deepcopy(self)
class Container(Ball):
    def __init__(self, radius=10, mass=100):
        super().__init__(); self._radius = radius; self._mass = mass; self._count = -1