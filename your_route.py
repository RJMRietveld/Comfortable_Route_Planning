"""
This file is part of the Bachelor End Project (BEP)

Comfortable  Route  Planning

from mechanical engineering of the Delft university of technology.

Authors: Bas de Jager, Maurits Neele, Mike Rietveld and Julian waas
"""

"""
Import modules
"""
from IPython import get_ipython
import route_generator_definitions
import generate_least_sickening_route_definitions

"""
Reset variable explorer
"""
get_ipython().magic('reset -sf')


"""
Input your route_name, origin_coordinates and destination_coordinates
"""
name_of_route = "#route_name"
origin_coordinates_input = 51.876533, 4.530575
destination_coordinates_input = 51.932910, 4.456871


"""
Structure of the algorithm

your_route
    -> route_generator_definitions
    -> generate_least_sickening_route_definitions
        -> route_generator_definitions
        -> route_calculations_definitions
        -> speed_and_acceleration_definitions
        -> motion_sickness_definitions
        -> LtiManip.py
"""

"""
Run definitions of this module
"""
G, G_4, routes = route_generator_definitions.run_route_generator_definitions(\
    origin_coordinates_input, destination_coordinates_input)

dataset_results, dataset_end = generate_least_sickening_route_definitions.\
    run_generate_least_sickening_route_definitions(\
    G, G_4, routes, name_of_route, origin_coordinates_input,\
    destination_coordinates_input)