import scipy.integrate as integrate
import numpy as np
from actin import actin
from myosin import myosin
from okt3 import okt3
from random import random
import random as rndm
from simulation_values import simulation_values
from mol_type import mol_type
from point_near import point_near
from norm import norm
import math

class synapse_simulation:
    """A class that runs the neccessary steps of an imune synapse
    simulation
    """

    def __init__(self,speeds):
        self.width = 100.
        self.height = 100.
        # self.mol_all is a list of molecule classes. All the
        # molecules in the simulation are contained in here.
        self.mol_all = list([])
        self.step_num = int(0)
        self.myosin_add = int(0)
        self.actin_add = int(0)
        # self.time is used in the numerical solving of the
        # differential equations that govern both population growth
        # and movement.
        self.time = np.array([])
        self.how_many = int(0)
        self.myo_count = int(0)
        self.act_count = int(0)
        self.radii_matrix = np.array([])
        self.distances_matrix = np.array([])
        if (simulation_values.ok_count > 1):
            self.band = int(0)
        else:
            self.band = "n/a"
        if any(speeds):
            simulation_values.term_velocity_act = speeds[0]
            simulation_values.term_velocity_myo = speeds[1]

    def step(self):
        """Each step of the simulation involves all the functions in
        self.step(). The simulation class calls this.
        """
        self.verhulst()  # Figures out how many molecules to add
        self.generate_molecules()
        self.move()
        self.age_molecules()
        self.set_valid()
        # Deletes invalid molecules
        self.mol_all = [i for i in self.mol_all if i.valid]
         
    def molecule_count(self):
        self.how_many = len(self.mol_all)

    def age_molecules(self):
        for molecule in self.mol_all:
            dtime = simulation_values.time_step
            if molecule.t == mol_type.MYOSIN:
                dtime = simulation_values.time_step*simulation_values.myosin_aging
            if simulation_values.overactive_actin:
                if molecule.t == mol_type.ACTIN:
                    dtime = simulation_values.time_step*simulation_values.actin_aging
            molecule.set_age(dtime)

    def time_array(self,resolution):
        self.time = np.linspace( \
            self.step_num*simulation_values.time_step, \
            (self.step_num+1)*simulation_values.time_step, \
            resolution)


    def initialize(self):
        """Initializes the molecules for the start of the synapse
        simulation. All myosin and actin are randomly places to
        simulate t=0 where the T-cell has not been introduced to OKT3.
        """

        if simulation_values.ok_count == 3:
            # Generates 3 OKT3
            a = okt3(25.,25.)
            b = okt3(75.,25.)
            c = okt3(50.,75.)
            self.mol_all.append(a)
            self.mol_all.append(b)
            self.mol_all.append(c)
        elif simulation_values.ok_count == 2:
            # Generates 2 OKT3
            for i in range(2):
                mol = okt3(self.width/4.*(2.*(i+.5)),self.height/2.)
                self.mol_all.append(mol)
        elif simulation_values.ok_count == 1:
            # Generates 1 OKT3
            a = okt3(50.,50.)
            self.mol_all.append(a)

        if simulation_values.ok_count == 1:
            # Generates actin radially for cortex
            for i in range(simulation_values.num_act):
                angle = random()*2.*np.pi
                mol = actin( \
                    ((simulation_values.cortex_rad+(random()-random())*4.+1.5))*np.cos(angle)+50., \
                    ((simulation_values.cortex_rad+(random()-random())*4.+1.5))*np.sin(angle)+50.)
                self.mol_all.append(mol)
        else:
            # Generates appropriate amount of myosin and actin for new
            # population dynamics
            for i in range(simulation_values.num_act):
                mol = actin(random()*self.width,random()*self.height)
                self.mol_all.append(mol)
        #for i in range(int(simulation_values.carrying_capacity_myo*(5./8.))):
        for i in range(simulation_values.num_myo):
            mol = myosin(random()*self.width,random()*self.height)
            self.mol_all.append(mol)
        
        self.molecule_count()

    def set_valid(self):
        """Sets the validity of each molecule by calling
        molecule.set_invalid() for each molecule. Validity determines
        if the molecule is deleted at the end of the step or
        not. Determined by position and age.
        """

        for molecule in self.mol_all:
            x = molecule.pos[0]
            y = molecule.pos[1]
            pos_val = ((0 < x < self.width) and (0 < y < self.height))

            age_val = ((molecule.age < simulation_values.decay_age + \
                        simulation_values.decay_age/3 * \
                        (random()*2-1.)) \
                       or (molecule.t == mol_type.OKT3))

            molecule.set_valid((pos_val) and (age_val))


    def generate_molecules(self):
        """Generates the number of molecules prescribed by
        self.verlhurst()
        """
        
        for i in range(self.myosin_add):
            self.mol_all.append(myosin(random()*self.width,random()*self.height))
        for i in range(self.actin_add):
            where = random()
            if (where <= simulation_values.near_okt3_odds):
                # Generates actin near OKT3 by choosing points near
                # the first ok_count molecules which are the OKT3s
                which = int(random()*float(simulation_values.ok_count))
                okt3 = self.mol_all[which]
                p = point_near(okt3.pos,5.)
                self.mol_all.append(actin(p[0],p[1]))
            else:
                # Generates actin randomly on field
                if simulation_values.ok_count == 1:
                    angle = random()*2.*np.pi
                    self.mol_all.append( \
                        actin((simulation_values.cortex_rad+3.)*random()*np.cos(angle)+50., \
                              (simulation_values.cortex_rad+3.)*random()*np.sin(angle)+50.) \
                        )
                else:      
                    self.mol_all.append(actin(random()*self.width,random()*self.height))
        self.molecule_count()
            

    def get_radii(self,positions):
        """This function sets self.radii_matrix and self.distances_matrix
        to the appropraite values to be used in the equation of mation
        """

        # The radius matrix's dimensions are
        # [self.how_many,self.how_many,2]. This is because each radius
        # is actually a vector of length 2. The distance matrix's
        # dimensions are [self.how_many,self.how_many] because it is
        # an array of scalars.
        rad_mat = np.zeros([self.how_many,self.how_many,2])
        dist_mat = np.zeros([self.how_many,self.how_many])
        for i in range(self.how_many):
            for j in range(self.how_many):
                # To avoid having to do excess calculations, if
                # the row number is larger than the column number,
                # the new radius and distance is set to the negative
                # of the one with the opposite index
                if (i > j):
                    rad_mat[i,j] = -rad_mat[j,i]
                    dist_mat[i,j] = dist_mat[j,i]
                elif (i == j):
                    rad_mat[i,j] = np.zeros([2])
                    dist_mat[i,j] = float(0)                
                elif (i != j):
                    rad_mat[i,j] = \
                                 positions[j] - \
                                 positions[i]
                    dist_mat[i,j] = norm(rad_mat[i,j])
                                  
        return rad_mat,dist_mat


    def move(self):
        """This method moves all the molecules in a manner backed by
        reasonably sound mathematics.
        """
        
        # Position "differential." This number is the length of the
        # vector of displacement for each molecule every time it
        # moves. Essentially, there is a terminal velocity
        # (simulation_values.term_velocity) and the rest of this
        # method figures out the unit vector that dictates it's direction.
        dv_act = simulation_values.term_velocity_act \
             / (simulation_values.num_steps*simulation_values.resolution_motion)
        dv_myo = simulation_values.term_velocity_myo \
             / (simulation_values.num_steps*simulation_values.resolution_motion)
        
        def mover():
            S_i = np.zeros(self.how_many*2)
            for i in range(self.how_many):
                S_i[2*i] = self.mol_all[i].pos[0]
                S_i[2*i+1] = self.mol_all[i].pos[1]
            eom(S_i)
            
        def eom(W):
            w = W.reshape([self.how_many,2])
            radii_matrix, distances_matrix = self.get_radii(w)
            dist_max2 = simulation_values.attract_upper
            dist_min2 = simulation_values.attract_lower
            
            for i in range(self.how_many):
                rx = 0.
                ry = 0.
                r0 = bool(0)
                if (self.mol_all[i].t == mol_type.OKT3):
                    r0 = bool(1)
                else:
                    for k in range(self.how_many):
                        if (i == k):
                            pass
                        elif (distances_matrix[i,k] > dist_max2):
                            pass
                        # For when there's only one OKT3
                        elif ((simulation_values.ok_count == 1) and \
                              (self.mol_all[i].t == mol_type.ACTIN) and \
                              (self.mol_all[k].t == mol_type.OKT3) and \
                              (distances_matrix[i,k] > simulation_values.cortex_rad)):
                           r0 = bool(1)
                           self.mol_all[i].age = self.mol_all[i].age*.1
                           break
                        elif (((self.mol_all[i].t == mol_type.ACTIN) and \
                               (self.mol_all[k].t == mol_type.ACTIN)) and \
                              (distances_matrix[i,k] < simulation_values.repulsion_radius_act)):
                            dm = distances_matrix[i,k] 
                            rx += -radii_matrix[i,k][0]/(dm/simulation_values.repuls_factor)
                            ry += -radii_matrix[i,k][1]/(dm/simulation_values.repuls_factor)
                        elif (((self.mol_all[i].t == mol_type.MYOSIN) and \
                               (self.mol_all[k].t == mol_type.MYOSIN)) and \
                              (distances_matrix[i,k] < simulation_values.repulsion_radius_myo)):
                            dm = distances_matrix[i,k]                            
                            rx += -radii_matrix[i,k][0]/(dm/simulation_values.repuls_factor)
                            ry += -radii_matrix[i,k][1]/(dm/simulation_values.repuls_factor)
                        elif (distances_matrix[i,k] < dist_min2):
                            pass
                        elif ((self.mol_all[i].t == mol_type.MYOSIN) and \
                              (self.mol_all[k].t == mol_type.ACTIN)):
                            dm = distances_matrix[i,k]
                            rx += radii_matrix[i,k][0]/dm
                            ry += radii_matrix[i,k][1]/dm
                        elif ((self.mol_all[i].t == mol_type.ACTIN) and \
                              (self.mol_all[k].t == mol_type.MYOSIN)):
                            dm = distances_matrix[i,k]
                            rx += radii_matrix[i,k][0]/dm
                            ry += radii_matrix[i,k][1]/dm
                    if ((rx != 0) or (ry != 0)):
                        rm = norm(np.array([rx,ry]))
                        if (rm > 1.):
                            rx = rx/math.sqrt(rm)
                            ry = ry/math.sqrt(rm)
                if r0:
                    vx = 0.
                    vy = 0.
                if (self.mol_all[i].t == mol_type.ACTIN):
                    vx = rx*dv_act
                    vy = ry*dv_act
                if (self.mol_all[i].t == mol_type.MYOSIN):
                    vx = rx*dv_myo
                    vy = ry*dv_myo                
                self.mol_all[i].set_position([self.mol_all[i].pos[0]+vx, \
                                              self.mol_all[i].pos[1]+vy])
                    
        for i in range(simulation_values.resolution_motion):
            mover()


    def verhulst(self):
        """Finds the appropriate number of myosin and actin to add to
        the simulation at every step. Does so by using two related verhulst
        differential equations.
        """
        
        # Declares instance variables for initial population
        # conditions
        Pm_i = 0.
        Pa_i = 0.
        
        for molecule in self.mol_all:
            if (molecule.t == mol_type.MYOSIN):
                Pm_i += 1.
            if (molecule.t == mol_type.ACTIN):
                Pa_i += 1.
        
        def population2(P,t):
            K_myo = simulation_values.carrying_capacity_myo
            K_act = simulation_values.carrying_capacity_act
            Pm = P[0]
            Pa = P[1]
            dPdt = np.array(\
                [simulation_values.myosin_genesis_rate*Pm*(1.-(Pm)/K_myo),\
                 simulation_values.actin_genesis_rate*Pa*(1.-(Pa)/K_act)])
            return dPdt

        self.time_array(simulation_values.resolution_population)
        # Initial conditions to pass to differential equation solver
        P_i = [Pm_i,Pa_i]
        P = integrate.odeint(population2,P_i,self.time)
        # Following two lines figure out how many more myosin and
        # actin need to be added
        self.myo_count = P[-1,0]
        self.act_count = P[-1,1]
        self.myosin_add = int(P[-1,0] - Pm_i)
        self.actin_add = int(P[-1,1] - Pa_i)


    def dist_to_band(self):
        if (simulation_values.ok_count > 1):
            self.band = 0
            for molecule in self.mol_all:
                self.band += math.fabs(self.height/2 - molecule.pos[1])
            self.band = self.band/self.how_many
        else:
            self.band = "n/a"
