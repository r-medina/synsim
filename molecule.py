import numpy as np
from random import random

class molecule:
    """Simple molecule class."""

    def __init__(self,position):
        self.pos = position
        self.valid = bool(1)
        self.age = 40.*random()
        self.vel = np.zeros(2)
        self.hist = self.pos
        
    def set_valid(self,yn):
        self.valid = yn

    def set_age(self,time_to_add):
        self.age += time_to_add

    def set_position(self,new_position):
        self.pos = new_position
        #self.save_hist(self.pos)

    def set_velocity(self,new_velocity):
        self.vel = new_velocity

    # Saves all of the positions of a given molecule over its lifetime
    def save_hist(self,position):
        self.hist = np.vstack([self.hist,position])
