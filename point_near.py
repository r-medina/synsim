from random import random
import math

def point_near(point,factor):
    return point[0]+(random()*2.*factor-factor),point[1]+random()*2.*factor-factor
