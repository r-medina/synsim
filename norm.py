import numpy as np

def norm(vector):
    """Finds the length of a vector"""
    if all(vector == np.zeros(len(vector))):
        mag = float(0)
    else:
        mag = float((np.dot(vector,vector)))
    return mag
