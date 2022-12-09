
class Event(tuple):
    """ EVENT CLASS |
        A tuple of 5 elements (ball_A, ball_B, count_A, count_B, dt).
        < PARAMETERS >
        -> ball_A (int): The first ball in impending collision.
        -> ball_B (int): The second ball in impending collision.
        -> count_A (int): The number of collisions the first ball
                encountered prior to this impending collision calculation.
        -> count_B (int): The number of collisions the second ball
                encountered prior to this impending collision calculation.
        -> dt (float): The global time this collision will happen on.
        RETURNS
            pair (list(int)): a list of the two balls or a ball and container
                              involved in the collision.
            ball_A (int): the index of first ball in impending collision.
            ball_B (int): the index of second ball in impending collision.
            count_A (int): the number of collisions the first ball encountered
                           prior to this impending collision calculation.
            count_B (int): the number of collisions the second ball encountered
                           prior to this impending collision calculation.
            dt (float): the global time(step) this collision will happen on.
    """
    def pair(self):         return [self[0], self[1]]   # Returns pair object.
    def ball_A(self):       return self[0]              # Returns Ball A.
    def ball_B(self):       return self[1]              # Returns Ball B.
    def count_A(self):      return self[2]              # Returns Ball A count.
    def count_B(self):      return self[3]              # Returns Ball B count.
    def dt(self):           return self[4]              # Returns minimum time.
    