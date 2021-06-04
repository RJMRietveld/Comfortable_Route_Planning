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
from scipy import integrate, signal
from LtiManip import ltimul
from scipy.signal import find_peaks
import warnings

"""
Run definitions of this module
"""
def run_motion_sickness_definitions(v_accurate, a_lat, a_long,\
    generated_interpolation_step):
    
    t, dt, v_t, a_long_t, a_lat_t = time_dependent(v_accurate, a_lat, a_long,\
        generated_interpolation_step)
    
    filtered_acceleration_long, filtered_acceleration_lat = filter_data(a_long_t,\
        a_lat_t, t)
    
    misc, misc_max, misc_peaks = filter_data_oman(filtered_acceleration_long,\
        filtered_acceleration_lat, t)
    
    msdv = calculate_msdv_integral(filtered_acceleration_lat,\
        filtered_acceleration_long, t)
    
    cwa, cwa_step = cumulative_weighted_accelerations(t, dt,\
        filtered_acceleration_lat, filtered_acceleration_long)
    
    return(msdv, cwa, cwa_step, t, misc, misc_max, misc_peaks, v_t, a_long_t,\
           a_lat_t, filtered_acceleration_long, filtered_acceleration_lat)

"""
Definitions of this module
"""
def time_dependent(v_x, a_lat_x, a_long_x, ds_array):
    """
    Function transforming arrays dependend on distance to arrays dependend on time.
    
    Parameters
    ----------
    v_x : 
        Speed array with data dependent on x.
    a_lat_x : 
        Lateral acceleration array with data dependent on x.
    a_long_x : 
        Longitudal acceleration array with data dependent on x.
    ds_array : 
        Array containing the distance between datapoints.
    
    Returns
    -------
    t :
        Time array.
    v_t :
        Speed array with data dependent on t.
    a_lat_t :
        Lateral acceleration array with data dependent on t.
    a_long_t
        Longitudal acceleration array with data dependent on t.
    """
    t_datapoint = np.zeros(len(v_x))
    
    for i in range(len(v_x)-1):
        if v_x[i] == 0:
            t_datapoint[i+1] = t_datapoint[i]
        elif 2*a_long_x[i]*ds_array[i]+(v_x[i])**2 < 0:
            t_datapoint[i+1] = t_datapoint[i] + ds_array[i]/(v_x[i])
        elif a_long_x[i] == 0:
            t_datapoint[i+1] = t_datapoint[i] + ds_array[i]/(v_x[i])

        else:
            t_datapoint[i+1] = t_datapoint[i] + (np.sqrt(2*a_long_x[i]*ds_array\
                [i]+(v_x[i])**2) - v_x[i])/a_long_x[i]

    dt = 0.01                                   
    t = np.linspace(0, t_datapoint[-1], round(t_datapoint[-1]/dt))
    
    v_t = np.zeros(len(t))
    a_long_t = np.zeros(len(t))
    a_lat_t = np.zeros(len(t)) 
    k = 0  
    for i in range(len(t)):   
        v_t[i] = v_x[k]
        a_long_t[i] = a_long_x[k]
        a_lat_t[i] = a_lat_x[k]
        if t[i] >= t_datapoint[k]:   
            k = k + 1
    
    return(t, dt, v_t, a_long_t, a_lat_t)
    
def filter_data(a_long_t, a_lat_t, t):
    """
    Motion_sickness_filter data using filter system.
    
    Parameters
    ----------
    filter_no : 
        Choose which filter you want to use from motion_sickness_filter function.
    data : array_like
        An input array describing the input at each time `T`
        (interpolation is assumed between given times).  If there are
        multiple inputs, then each column of the rank-2 array
        represents an input.  If U = 0 or None, a zero input is used.
    t : array_like
        The time steps at which the input is defined and at which the
        output is desired.  Must be nonnegative, increasing, and equally spaced.
    
    Returns
    -------
    filtered_data :
        Array with filtered data.
    """
    filter_no = 8
    filter_no = filter_no - 1
    
    H_transfer = motion_sickness_filter()[filter_no]
    n,d = H_transfer.num, H_transfer.den 
    
    filtered_acceleration_long = signal.lsim((n,d), a_long_t, t)[1]
    filtered_acceleration_lat = signal.lsim((n,d), a_lat_t, t)[1]
    
    return(filtered_acceleration_long, filtered_acceleration_lat)
  
def filter_data_oman(a_long_t, a_lat_t, t):
    """
    Motion_sickness_filters data using filter system.
    
    Parameters
    ----------
    filter_no : 
        choose which filter you want to use from motion_sickness_filter function
    data : array_like
        An input array describing the input at each time `T`
        (interpolation is assumed between given times).  If there are
        multiple inputs, then each column of the rank-2 array
        represents an input.  If U = 0 or None, a zero input is used.
    t : array_like
        The time steps at which the input is defined and at which the
        output is desired.  Must be nonnegative, increasing, and equally spaced.
    
    Returns
    -------
    filtered_data :
        Array with filtered data.
    """
    conflict = np.sqrt(a_long_t**2+a_lat_t**2) #+a_lat_t**2
   
    Theta_1 = 526.462778190096
    Theta_2 = 57.1077822064923
    Theta_3 = 6.90548977423296 
    
    N_SP = [Theta_3]
    D_SP = [Theta_1**2, 2*Theta_1,1]
    N_FP= [1]
    D_FP = [Theta_2**2,2*Theta_2,1]
    
    slowpathout= signal.lsim((N_SP,D_SP),conflict,t)[1] 
    inputfastpath = np.multiply(conflict,slowpathout) 
    misc = signal.lsim((N_FP,D_FP),inputfastpath,t)[1] 
    misc_max = max(misc)
    misc_peaks, _ = find_peaks(misc, height=6)

    return(misc, misc_max, misc_peaks)

def motion_sickness_filter():

    warnings.simplefilter("ignore")        
    
    """
    Notes
    ------
    Code provided by Tugrul Irmak
    """
    dt = 0.01
    e = 10**(-10)
    inf= 10**10 
    pi =  np.pi
    cnst = np.array([[2*pi*0.4, 1/(2**0.5), 2*pi*100, 1/(2**0.5), 2*pi*8, 2*pi*8,\
        0.63, 0, e, 0, e],
    [2*pi*0.4, 1/(2**0.5), 2*pi*100, 1/(2**0.5), 2*pi*2, 2*pi*2, 0.63, 0, e, 0, e],
    [2*pi*0.4, 1/(2**0.5), 2*pi*100, 1/(2**0.5), 2*pi*1, 2*pi*1, 0.63, 0, e, 0, e],
    [2*pi*0.08, 1/(2**0.5), 2*pi*0.63, 1/(2**0.5), 2*pi*(1/(2*dt)), 2*pi*0.25,\
         0.86,2*pi*0.06, 0.8,  2*pi*0.1, 0.8],
   [2*pi*0.4, 1/(2**0.5), 2*pi*100, 1/(2**0.5), e, 0, e,  2*pi*3.75, 0.91,\
        2*pi*5.3, 0.91],
    [2*pi*0.4, 1/(2**0.5), 2*pi*100, 1/(2**0.5), 2*pi*12.5, 2*pi*12.5, 0.63,\
         2*pi*2.37, 0.91,  2*pi*3.3, 0.91],
    [2*pi*0.7943, 1/(2**0.5), 2*pi*100, 1/(2**0.5), 2*pi*5.68, 2*pi*5.684, 0.5,\
         0, e, 0, e],
    [2*pi*0.02,  1/(2**0.5), 2*pi*0.63, 1/(2**0.5), inf, 2*pi*0.25, 0.86, inf,\
         1, inf, 1]]) 

    motion_sickness_filters = list([] for x in range(8)) 
    
    # This produces filters in the followign form:
    # motion_sickness_filter{1} - Fore-aft backrest vibration ISO 2631-1
    # motion_sickness_filter{2} - Fore-aft and lateral Seat vibration ISO 2631-1
    # motion_sickness_filter{3} - Roll, pitch, yaw Rotational seat vibration 
    #   ISO 2631-1
    # motion_sickness_filter{4} - Vertical Motion sickness ISO 2631-1
    # motion_sickness_filter{5} - Vertical Head vibration ISO 2631-1
    # motion_sickness_filter{6} - Vertical Seat vibration ISO 2631-1
    # motion_sickness_filter{7} - Tall Building vibration ISO 2631-2
    # motion_sickness_filter{8} - Griffin lateral vibrations 
    #Data is from http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.324.
    #   3870&rep=rep1&type=pdf
    #"Design of Digital motion_sickness_filters for Frequency Weightings 
    #   Required for Risk Assessments of Workers Exposed to Vibration"
    
    for x in range(len(motion_sickness_filters)): 
        Hl = ltimul(cnst[x][2]**2,[1,cnst[x][2]/cnst[x][3], cnst[x][2]**2])
        Hh = ltimul([1,0,0],[1,cnst[x][0]/cnst[x][1], cnst[x][0]**2 ])
        Ht = ltimul([(cnst[x][5]**2)/cnst[x][4],(cnst[x][5]**2)],[1,cnst[x][5]\
            /cnst[x][6], cnst[x][5]**2])
        Hs = ltimul([1,(cnst[x][7]/cnst[x][8]),cnst[x][7]**2],[1, (cnst[x][9]\
            /cnst[x][10]), cnst[x][9]**2])
    
        if x ==3 or x == 5:
            motion_sickness_filters[x] = Hh*Hl*Ht*Hs
        if x <= 2:
            motion_sickness_filters[x] = Hh*Hl*Ht
        if x == 6:
            motion_sickness_filters[x] = Hh*Hl*Ht
        if x == 7:
            motion_sickness_filters[x] = 1.72*0.55*Hh*Hl*Ht        
        else:
            motion_sickness_filters[x] = Hh*Hl*Hs
    
    motion_sickness_filters = motion_sickness_filters
    
    warnings.simplefilter("default")        
    
    return(motion_sickness_filters)

"""
Definitions from below define two other methods to indicate motion sickness,
the Motion_Sickness_Dosage_Value (MSDV) and Cumulative_Weighted_Accelerations
(CWA) method.
"""

def calculate_msdv_integral(filtered_acceleration_lat,filtered_acceleration_long,\
    t):
    """
    Function calculating MSDV.
    
    Htike, Z., Papaioannou, G., Siampis, E., Velenis, E., & Longo, S. (2020, June). 
    Motion Sickness Minimisation in Autonomous Vehicles Using Optimal Control. 
    In International Conference on Robotics in Alpe-Adria Danube Region (pp. 
    275-282). Springer, Cham.
    
    Parameters
    ----------
    t : 
        Time array.
    a_lat_filtered : 
        Array with filtered lateral accelerations.
    a_long_filtered : 
        Array with filtered longitudal accelerations.
    
    Returns
    -------
    MSDV : 
        Returns msdv value.
    
    Notes
    ------
    MSDV not that suitable for comparing routes.
    Sse Cumulative_Weighted_Accelerations instead.
    """
    integral_lat = integrate.simps((filtered_acceleration_lat)**2, t)
    integral_long = integrate.simps((filtered_acceleration_long)**2, t)
    
    msdv = np.sqrt(integral_lat) + np.sqrt(integral_long)
    
    return(msdv)

def cumulative_weighted_accelerations(t, dt, filtered_acceleration_lat,\
    filtered_acceleration_long):
    """
    Function calculating cumulative weighted accelerations.
    
    Parameters
    ----------
    a_lat_filtered : 
        Array with weighted lateral acceleration.
    a_long_filtered : 
        Array with longitudal accelerations.
    t : 
        Time array.
    
    Returns
    -------
    cwa : 
        Cumulative weighted accelerations.
    """
    
    cwa = np.zeros(len(t))
    cwa_step = np.zeros(len(t))
    
    for i in range(len(t)-1):
        cwa[i+1] = cwa[i] + ((filtered_acceleration_lat[i])**2 + (\
            filtered_acceleration_long[i])**2)*dt
        cwa_step[i+1] = cwa[i+1]-cwa[i]
    
    return(cwa,cwa_step)