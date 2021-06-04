"""
This file is part of the Bachelor End Project (BEP)

Comfortable  Route  Planning

from mechanical engineering of the Delft university of technology.

Authors: Bas de Jager, Maurits Neele, Mike Rietveld and Julian waas
"""

"""
Import modules
"""
import numpy as np

"""
Run definitions of this module
"""
def run_speed_and_acceleration_definitions(generated_interpolation_step,\
    K, max_speed_per_node):
    
    v_accurate, a_lat, a_long = accurate_speed_profile(\
        generated_interpolation_step, K, max_speed_per_node)
    
    return(v_accurate, a_lat, a_long)

"""
Definitions of this module
"""
def accurate_speed_profile(generated_interpolation_step, K, max_speed_per_node):
    """
    Function creating accurate speed profile.
    
    Parameters
    ----------
    max_speed_per_node : 
        Speed array of legal speed limit at each datapoint.
    K : 
        Curvature at each datapoint.
    generated_interpolation_step : 
        Distance between datapoints.
    
    Returns
    -------
    v_accurate : 
        Speed array.
    a_long : 
        Longitudal accelerations array.
    a_Lat : 
        Lateral accerations array.
    """
   
    max_speed_recalculated, a_lat = max_speed(max_speed_per_node, K)
    
    v_accurate, a_long = longitudal_constant_acceleration(max_speed_recalculated,\
        generated_interpolation_step)
    
    return(v_accurate, a_lat, a_long)

def longitudal_constant_acceleration(max_speed_recalculated,\
    generated_interpolation_step): 
    """
    Function returning speed profile, calculated with constant accelerations
    
    Parameters
    ----------
    max_speed_recalculated : 
        Speed array with max speed based on constant longitudinal acceleration.
    generated_interpolation_step : 
        Distance between datapoints.
   
    Returns
    -------
    v_accurate : 
        Accurate speed array with accelerations.
    a_long : 
        Longitudinal accelerations at each point.
    """
    
    a = 1.5   
    da = -1.5   
    
    v = max_speed_recalculated
    
    a_l = np.zeros(len(v))
        
    for i in range(len(v)-1):
        if  (v[i+1]-.5) <= v[i] <= (v[i+1]+0.5):   
            v[i+1] = v[i]
        elif v[i] <= v[i+1]:                             
            v[i+1] = np.sqrt((v[i])**2 + 2*a*generated_interpolation_step[i])
            a_l[i] = a
    
    for i in range(len(v)-1,0,-1):
        if  (v[i-1]-0.1) <= v[i] <= (v[i-1]+0.1):   
                v[i-1] = v[i]
        elif v[i-1] >= v[i]:
                v[i-1] = np.sqrt((v[i])**2 - 2*da*generated_interpolation_step[i])
                a_l[i] = da
    
    v_accurate = v
    a_long = a_l
    
    return(v_accurate,a_long)

def max_cornering_speed(K_i):  
    """
    Function calculating maximum cornering speed dictated by friction and radius
    Can be used to compare the maximum speed dictated by turns and legal speed limits.
    
    Parameters
    ----------
    K_i : 
        Local radius turn.
    
    Returns
    -------
    v_max_turn :
        Maximum speed in turn.
    
    Notes:
        Kapania, Nitin & Subosits, John & Gerdes, J.. (2016). A Sequential Two-
        Step Algorithm for Fast Generation of Vehicle Racing Trajectories. 
        Journal of Dynamic Systems, Measurement, and Control. 
        138. 10.1115/1.4033311. 
        
    Constraint by a certain acceleration of 1.5 m/s^2
    """
    
    a_lat_max = 1.5   
     
    if abs(K_i) > 0:                                               
        v_max_turn = np.sqrt(a_lat_max*(1/abs(K_i)))      
    else:
        v_max_turn = 80           
    
    return(v_max_turn)

def max_speed(max_speed_per_node, K): 
    """
    Function comparing max legal speed with max speed through turns and returning 
        minimum value
    
    Parameters
    ----------
    max_speed_per_node : 
        Array containing max node speegenerated_interpolation_step.
    K : 
        Array containing Curvature at each point.
    
    Returns
    -------
    max_speed_recalculated: 
        Max speed on each point.
    a_lat : 
        Lateral accelerations at each point.
    """
    V_turn_max = np.zeros(len(max_speed_per_node))
    a_lat = np.zeros(len(max_speed_per_node))
    for i in range(len(max_speed_per_node)):
        V_turn_max[i] = max_cornering_speed(K[i])
        
        if V_turn_max[i] <= max_speed_per_node[i]:
                max_speed_per_node[i]=V_turn_max[i]
                a_lat[i] = K[i]*(V_turn_max[i])**2 
        else:
            a_lat[i] = K[i]*(max_speed_per_node[i])**2      
    
    max_speed_recalculated = max_speed_per_node
    
    max_speed_recalculated[0] = 0
    max_speed_recalculated[-1] = 0
    
    return(max_speed_recalculated, a_lat)