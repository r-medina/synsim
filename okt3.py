import numpy as np
from molecule import molecule
from mol_type import mol_type

class okt3(molecule):
    def __init__(self,x,y):
        molecule.__init__(self,np.array([x,y]))
        self.t = mol_type.OKT3
        
