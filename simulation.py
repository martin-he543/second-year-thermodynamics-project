import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import itertools as it
import time
import random

import ball as bl
import event as ev

class Simulation():
    """
    CLASS: Simulation
    Simulates rigid collisions in 2D within a circular container.

    Attributes:  MUST CHANGE
        N_ball (int, optional): Number of balls in the container.
        r_container (float, optional): Radius of container.
        r_ball (float, optional): Radius of balls in the container (if all
            balls have the same radius)
        m_ball (float, optional): Mass of balls in the container (if all balls
            have the same mass)
        random_pos (boolean, optional): If True, all balls are initialised with
            random positions.
        random_speed_range (numpy.ndarray of float, optional): Gives the range 
            of x- and y- speeds to be chosen randomly from a uniform 
            distribution [-random_speed_range, random_speed_range]
    """
    
    def __init__(
        self,
        ball_N=1,
        ball_r=1,
        ball_m=1,
        container_r=10,
        ball_random_pos=True,
        ball_random_speed_range=True,
    ):
        self._container = 10
        self._ball = []
        
    def __repr__(self):
        """
        Gives the representation format for simulation class in console.
        """
        return
    
    def __str__(self):
        """
        Gives the string format for simulation class in print statements.
        """
        return
    
    
    ## PROPERTIES OF THE SYSTEM
##.................................................................................................................
##.PPPPPP.....PRRRRRR..........OOO......OPPPPPP.....PEEEEEEEEE.ERRRRR.....RRTTTTTTTTTT.TII..IEEEEEEEE.....SSSS.....
##.PPPPPPPPP..PRRRRRRRRR....OOOOOOOOO...OPPPPPPPP...PEEEEEEEEE.ERRRRRRRR..RRTTTTTTTTTT.TII..IEEEEEEEE...SSSSSSSS...
##.PPPPPPPPPP.PRRRRRRRRRR..OOOOOOOOOOO..OPPPPPPPPP..PEEEEEEEEE.ERRRRRRRRR.RRTTTTTTTTTT.TII..IEEEEEEEE..ESSSSSSSSS..
##.PPPPPPPPPP.PRRRRRRRRRR..OOOOOOOOOOO..OPPPPPPPPPP.PEEEEEEEEE.ERRRRRRRRR.RRTTTTTTTTTT.TII..IEEEEEEEE..ESSSSSSSSS..
##.PPP...PPPP.PRRR...RRRR.ROOOO...OOOOO.OPPP...PPPP.PEEE.......ERRR...RRRR....TTTT.....TII..IEE........ESSS...SSS..
##.PPP....PPP.PRRR...RRRR.ROOO.....OOOO.OPPP...PPPP.PEEE.......ERRR...RRRR....TTTT.....TII..IEE........ESSSS.......
##.PPP...PPPP.PRRR..RRRRR.ROOO.....OOOO.OPPP...PPPP.PEEEEEEEE..ERRR..RRRRR....TTTT.....TII..IEEEEEEEE..ESSSSSSS....
##.PPPPPPPPPP.PRRRRRRRRRR.ROOO.....OOOO.OPPPPPPPPPP.PEEEEEEEE..ERRRRRRRRR.....TTTT.....TII..IEEEEEEEE...SSSSSSSSS..
##.PPPPPPPPPP.PRRRRRRRRR..ROOO.....OOOO.OPPPPPPPPP..PEEEEEEEE..ERRRRRRRRR.....TTTT.....TII..IEEEEEEEE.....SSSSSSS..
##.PPPPPPPP...PRRRRRRRR...ROOO.....OOOO.OPPPPPPPP...PEEE.......ERRRRRRRR......TTTT.....TII..IEE........ESS...SSSS..
##.PPP........PRRR..RRRR..ROOOO...OOOOO.OPPP........PEEE.......ERRR.RRRR......TTTT.....TII..IEE........ESSS...SSS..
##.PPP........PRRR..RRRRR..OOOOOOOOOOO..OPPP........PEEEEEEEEE.ERRR..RRRR.....TTTT.....TII..IEEEEEEEE..ESSSSSSSSS..
##.PPP........PRRR...RRRR..OOOOOOOOOOO..OPPP........PEEEEEEEEE.ERRR..RRRRR....TTTT.....TII..IEEEEEEEE..ESSSSSSSSS..
##.PPP........PRRR...RRRRR...OOOOOOO....OPPP........PEEEEEEEEE.ERRR...RRRR....TTTT.....TII..IEEEEEEEE...SSSSSSSS...
##.............................OOO........................................................................SSSS.....
##.................................................................................................................

    def N_ball(self):
        """
        METHOD: N_ball
        Returns:
            (int): the number of balls in the container.
        """
        return
    
    def pressure(self):
        """
        METHOD: pressure
        Gives pressure values for every __ container collisions.

        Returns:
        """
        return
    
    def temperature(self):
        """
        METHOD: temperature
        Gives temperature at all collision times.
        
        Returns:
        """
        return
    
    def average_temperature(self):
        """
        METHOD: average_temperature
        Gives the average temperature of 
        
        Returns:
        """
        return
    
    def kinetic_energy(self):
        """
        METHOD: kinetic_energy
        Gives total kinetic energy of the system for all collisions.
        
        Returns:
        """
        return
        
    def speed(self):
        """
        METHOD: speed
        Gives all speeds for all balls for all collisions.
        
        Returns:
        """
        return
    
    
    ## METHODS OF THE SYSTEM