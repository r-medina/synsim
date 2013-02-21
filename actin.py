import numpy as np
from molecule import molecule
from mol_type import mol_type

class actin(molecule):
    def __init__(self,x,y):
        molecule.__init__(self,np.array([x,y]))
        self.t = mol_type.ACTIN
