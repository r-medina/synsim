import random as rndm

#rndm.seed("ba")
rndm.seed("ab")

class simulation_values():
    """The various values that are used throughout the simulation"""

    num_steps = 100
    carrying_capacity = 400
    ok_count = 2
    percent_myosin = .43
    carrying_capacity_myo = carrying_capacity*percent_myosin
    carrying_capacity_act = carrying_capacity*(1-percent_myosin)
    num_myo = int(round(0.5*carrying_capacity_myo))
    num_act = int(round(0.1*carrying_capacity_act))
    overactive_actin = False
    overactive_myosin = False
    actin_genesis_rate = .14 #.07
    myosin_genesis_rate = .2 #.09
    total_time = 250.
    time_step = float(total_time)/float(num_steps)
    decay_age = 35.1
    myosin_aging = .7 #.933 # Myosin stays around longere
    actin_aging = .6 # For overactive actin
    resolution_population = 25
    resolution_motion = 8
    term_velocity_act = 120.  # Changeable
    term_velocity_myo = 200.  # Changeable
    near_okt3_odds = .25
    repulsion_radius_act = 9.
    repulsion_radius_myo = 4.
    repuls_factor = 9.
    attract_upper = 50.**2.
    attract_lower = 2.**2. 
    cortex_rad = 40.
