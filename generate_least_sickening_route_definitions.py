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
import matplotlib.pyplot as plt
import pandas as pd
import osmnx as ox
from datetime import datetime
import os
import route_calculations_definitions
import speed_and_acceleration_definitions
import motion_sickness_definitions

"""
Run definitions of this module
"""
def run_generate_least_sickening_route_definitions(G, G_4, routes, name_of_route,\
    origin_coordinates_input, destination_coordinates_input):
    
    dataset_end = get_data_per_route(G, G_4, routes)
    
    dataset_results = get_least_sickening_route(G, G_4, routes, dataset_end,\
        name_of_route, origin_coordinates_input, destination_coordinates_input)
    
    return(dataset_results) 

"""
Definitions of this module
"""
def get_data_per_route(G, G_4, routes):
    """
    Create the desired number of different routes, get the data of these routes,
    calculate the variables that are wanted and save all those variables in a 
    dataframe. 

    Parameters
    -------
    G : networkx.MultiDiGraph
        Input graph.
    G_4 :  networkx.MultiDiGraph
        Input graph without any transformations.
    routes : list 
        List of routes.
    
    Returns
    -------
    dataset_end : DataFrame
        A pandas dataframe containing all information of the different routes.
        It contains the time, cwa_total, cwa_mean, cwa_max, distance, total_time
        msdv, route, misc, misc_mean, misc_max and misc_peaks. 
    """
    time = []
    distance_compare = []
    total_time_compare = []
    final_route_compare = []
    plot_index = 0
    
    misc_compare = []
    misc_mean_compare = []
    misc_max_compare = []
    misc_peaks_compare = []
    misc_minus_each_other = []
    
    msdv_compare = []
    cwa_compare = []
    cwa_mean_compare = []
    cwa_max_compare = []
    
    for route in routes:   
        print("")
        print("Route: ", plot_index)
        
        generated_interpolation_step, K, max_speed_per_node, dataset, ss =\
            route_calculations_definitions.run_route_calculations_definitions(G, G_4, route)
            
        v_accurate, a_lat, a_long = speed_and_acceleration_definitions.\
            run_speed_and_acceleration_definitions(generated_interpolation_step,\
            K, max_speed_per_node)
        
        msdv, cwa, cwa_step, t, misc, misc_max, misc_peaks, v_t, a_long_t,\
            a_lat_t, filtered_acceleration_long, filtered_acceleration_lat =\
            motion_sickness_definitions.run_motion_sickness_definitions(\
            v_accurate, a_lat, a_long, generated_interpolation_step)
        
        time.append(t)
        distance_compare.append(ss[-1])
        total_time_compare.append(t[-1])
        final_route_compare.append(route)
        
        misc_compare.append(misc)
        misc_mean_route = np.mean(misc)
        misc_mean_compare.append(misc_mean_route)
        misc_max_compare.append(misc_max)
        misc_peaks_compare.append(misc_peaks)
        
        misc_minus_each =  misc_max - misc_mean_route
        misc_minus_each_other.append(misc_minus_each)
        
        msdv_compare.append(msdv)
        cwa_compare.append(cwa[-1])
        cwa_mean_compare.append(np.mean(cwa_step))
        cwa_max_compare.append(max(cwa_step))

        plot_index += 1
                
    data_end = {"time":time ,"cwa_total":cwa_compare,"cwa_mean":cwa_mean_compare,\
                "cwa_max":cwa_max_compare,"distance":distance_compare,"total_time":\
                total_time_compare,"msdv":msdv_compare,"route":final_route_compare,\
                "misc":misc_compare,"misc_mean":misc_mean_compare,"misc_max":\
                misc_max_compare,"misc_peaks":misc_peaks_compare}
    dataset_end_not_unique = pd.DataFrame(data_end,columns=['time', 'cwa_total',\
        'cwa_mean','cwa_max','distance','total_time','msdv','route','misc',\
        'misc_mean','misc_max','misc_peaks'])
    dataset_end = dataset_end_not_unique.drop_duplicates(subset = ["cwa_total"])

    return(dataset_end)

def get_least_sickening_route(G, G_4, routes, dataset_end, name_of_route,\
    origin_coordinates_input, destination_coordinates_input):
    """
    Determines the final set of routes, the fastest route, least sickening 
    route (based on the lowest MISC-mean) and the route with the lowest max(MISC).
    Plots and a dataset containing all the information are create. Furthermore,
    a folder will be created containing the plots and a text-file containing
    the information. 

    Parameters
    -------
    G : networkx.MultiDiGraph
        Input graph.
    G_4 :  networkx.MultiDiGraph
        Input graph without any transformations.
    routes : list 
        List of routes.
    dataset_end : DataFrame
        A pandas dataframe containing all information of the different routes.
        It contains the time, cwa_total, cwa_mean, cwa_max, distance, total_time
        msdv, route, misc, misc_mean, misc_max and misc_peaks. 
    name_of_route : str
        Name of the route.
        Latitudal and longitudal coordinates of the origin.
    destination_coordinates_input : tuple
        Latitudal and longitudal coordinates of the destination.
    
    Returns
    -------
    dataset_results : DataFrame
        A pandas dataframe containing all information of the three final routes.
        The fastest route, least sickening route (based on the lowest MISC-mean)
        and the route with the lowest max(MISC).
    dataset_end : DataFrame
        Unchanged pandas dataframe containing all information of the different
        routes. It contains the time, cwa_total, cwa_mean, cwa_max, distance, 
        total_time, msdv, route, misc, misc_mean, misc_max and misc_peaks. 
    """
    local = datetime.now()
    local_string = local.strftime("%m-%d-%Y_%Hh%Mm")
    
    dirName = '{}_{}'.format(local_string,name_of_route)
    
    try:
        os.mkdir(dirName)
    except FileExistsError:
        print("")
        print("Directory " , dirName ,  " already exists")

    real_shortest_route = dataset_end['total_time'].min()
    index_real_shortest_route = dataset_end.loc[dataset_end['total_time'] ==\
        real_shortest_route].index.item()
   
    #Shortest_route
    shortest_route = dataset_end.at[index_real_shortest_route,'route']
    shortest_route_time = dataset_end.at[index_real_shortest_route,'time']
    shortest_route_total_time = dataset_end.at[index_real_shortest_route,\
        'total_time']
    shortest_route_distance = dataset_end.at[index_real_shortest_route,'distance']
    shortest_route_misc = dataset_end.at[index_real_shortest_route,'misc']
    shortest_route_misc_mean = dataset_end.at[index_real_shortest_route,'misc_mean']
    shortest_route_misc_max = dataset_end.at[index_real_shortest_route,'misc_max']
 
    figure_shortest_route, axis_shortest_route = ox.plot_graph_route(G,\
        shortest_route, show=False, close=False)
    plt.savefig('{}/{}_{}_fastest route no title.png'.format(dirName,\
        local_string,name_of_route))
    figure_shortest_route.suptitle('Fastest route', fontsize=20, color="black")
    plt.savefig('{}/{}_{}_fastest route.png'.format(dirName,\
        local_string,name_of_route))
    plt.show()
    
    plt.figure()
    plt.title("MISC - fastest route")
    plt.xlabel("time (s)")
    plt.ylabel("MISC (-)")
    plt.plot(shortest_route_time, shortest_route_misc, 'b', label='MISC')
    plt.hlines(shortest_route_misc_mean, shortest_route_time[0],\
        shortest_route_time[-1], 'r', '--', label='MISC mean = {}'.format(round(\
        shortest_route_misc_mean,1)))
    plt.hlines(shortest_route_misc_max, shortest_route_time[0],\
        shortest_route_time[-1], 'g', ':', label='MISC max = {}'.format(round(\
        shortest_route_misc_max,1)))
    plt.legend()
    plt.grid()
    plt.savefig('{}/{}_{}_misc_fastest route.png'.format(dirName,local_string,\
        name_of_route))
    
    
    #Lowest_misc_max
    minimum_misc_max = dataset_end['misc_max'].min()
    index_minimum_misc_max_route = dataset_end.loc[dataset_end['misc_max'] ==\
        minimum_misc_max].index.item()
    
    minimum_misc_max_route = dataset_end.at[index_minimum_misc_max_route,'route']
    minimum_misc_max_route_time = dataset_end.at[index_minimum_misc_max_route,\
        'time']
    minimum_misc_max_route_total_time = dataset_end.at[\
        index_minimum_misc_max_route,'total_time']
    minimum_misc_max_route_distance = dataset_end.at[index_minimum_misc_max_route,\
        'distance']
    minimum_misc_max_route_misc = dataset_end.at[index_minimum_misc_max_route,\
        'misc']
    minimum_misc_max_route_misc_mean = dataset_end.at[\
        index_minimum_misc_max_route,'misc_mean']
    minimum_misc_max_route_misc_max = dataset_end.at[\
        index_minimum_misc_max_route,'misc_max']
 
    figure_minimum_misc_max_route, axis_minimum_misc_max_route = ox.\
        plot_graph_route(G, minimum_misc_max_route, show=False, close=False)
    plt.savefig('{}/{}_{}_minimum max(misc) route no title.png'.format(\
        dirName,local_string,name_of_route))
    figure_minimum_misc_max_route.suptitle('Minimum max(misc) route',\
        fontsize=20, color="black")
    plt.savefig('{}/{}_{}_minimum max(misc) route.png'.format(dirName,\
        local_string,name_of_route))
    plt.show()
    
    plt.figure()
    plt.title("MISC - minimum max(misc) route")
    plt.xlabel("time (s)")
    plt.ylabel("MISC (-)")
    plt.plot(minimum_misc_max_route_time, minimum_misc_max_route_misc, 'b',\
        label='MISC')
    plt.hlines(minimum_misc_max_route_misc_mean, minimum_misc_max_route_time[0],\
        minimum_misc_max_route_time[-1], 'r', '--', label='MISC mean = {}'.\
        format(round(minimum_misc_max_route_misc_mean,1)))
    plt.hlines(minimum_misc_max_route_misc_max, minimum_misc_max_route_time[0],\
        minimum_misc_max_route_time[-1], 'g', ':', label='MISC max = {}'.format(\
        round(minimum_misc_max_route_misc_max,1)))
    plt.legend()
    plt.grid()
    plt.savefig('{}/{}_{}_misc_minimum max(misc) route.png'.format(dirName,\
        local_string,name_of_route))
    
    
    #Least_sickening_route
    minimum_misc_mean = dataset_end['misc_mean'].min()
    index_least_sickening_route = dataset_end.loc[dataset_end['misc_mean'] ==\
        minimum_misc_mean].index.item()
    
    least_sickening_route = dataset_end.at[index_least_sickening_route,'route']
    least_sickening_route_time = dataset_end.at[index_least_sickening_route,'time']
    least_sickening_route_total_time = dataset_end.at[index_least_sickening_route,\
        'total_time']
    least_sickening_route_distance = dataset_end.at[index_least_sickening_route,\
        'distance']
    least_sickening_route_misc = dataset_end.at[index_least_sickening_route,'misc']
    least_sickening_route_misc_mean = dataset_end.at[index_least_sickening_route,\
        'misc_mean']
    least_sickening_route_misc_max = dataset_end.at[index_least_sickening_route,\
        'misc_max']
 
    figure_least_sickening_route, axis_least_sickening_route =\
        ox.plot_graph_route(G, least_sickening_route, show=False, close=False)
    plt.savefig('{}/{}_{}_least sickening route no title.png'.format(dirName,\
        local_string,name_of_route))
    figure_least_sickening_route.suptitle('Least sickening route', fontsize=20,\
        color="black")
    plt.savefig('{}/{}_{}_least sickening route.png'.format(dirName,\
        local_string,name_of_route))
    plt.show()
    
    plt.figure()
    plt.title("MISC - least sickening route")
    plt.xlabel("time (s)")
    plt.ylabel("MISC (-)")
    plt.plot(least_sickening_route_time, least_sickening_route_misc, 'b',\
        label='MISC')
    plt.hlines(least_sickening_route_misc_mean, least_sickening_route_time[0],\
        least_sickening_route_time[-1], 'r', '--', label='MISC mean = {}'.\
        format(round(least_sickening_route_misc_mean,1)))
    plt.hlines(least_sickening_route_misc_max, least_sickening_route_time[0],\
        least_sickening_route_time[-1], 'g', ':', label='MISC max = {}'.format(\
        round(least_sickening_route_misc_max,1)))
    plt.legend()
    plt.grid()
    plt.savefig('{}/{}_{}_misc_least sickening route.png'.format(dirName,\
        local_string,name_of_route))
    
    #Dataframe
    result_route_number = [index_real_shortest_route, index_least_sickening_route,\
        index_minimum_misc_max_route]
    result_route_type = ["shortest","least sickening","minimum max(misc)"]
    result_time = [shortest_route_total_time, least_sickening_route_total_time,\
        minimum_misc_max_route_total_time]
    result_distance = [shortest_route_distance, least_sickening_route_distance,\
        minimum_misc_max_route_distance]
    result_misc_mean = [shortest_route_misc_mean, least_sickening_route_misc_mean,\
        minimum_misc_max_route_misc_mean]
    result_misc_max = [shortest_route_misc_max, least_sickening_route_misc_max,\
        minimum_misc_max_route_misc_max]    
    
    print("")
    print("MISC levels")
    symptom = ["No problems",
               "Slight discomfort but no specific symptoms",
               "Vague - Dizziness, warm, headache, stomach awareness, sweating,\
                  etc.",
               "Some - Dizziness, warm, headache, stomach awareness, sweating,\
                   etc.",
               "Medium - Dizziness, warm, headache, stomach awareness, sweating,\
                   etc.",
               "Severe - Dizziness, warm, headache, stomach awareness, sweating,\
                   etc.",
               "Some - Nausea",
               "Medium - Nausea",
               "Severe - Nausea", 
               "Retching - Nausea", 
               "Vomiting"]
    misc = [0,1,2,3,4,5,6,7,8,9,10]
    
    data_misc_values = {"symptom":symptom,"MISC":misc}
    dataset_misc_values = pd.DataFrame(data_misc_values,columns=['symptom',\
        'MISC'])
    print(dataset_misc_values)
    
    print("")
    data_results = {"number": result_route_number,"type":result_route_type,\
        "time (s)":result_time, "distance (m)":result_distance, "MISC mean (-)":\
        result_misc_mean, "MISC max (-)":result_misc_max}
    dataset_results = pd.DataFrame(data_results,columns=['number', 'type',\
        'time (s)', 'distance (m)', 'MISC mean (-)', 'MISC max (-)'])
    print(dataset_results)
    
    new_file = open('{}/{}_{}.txt'.format(dirName,local_string,name_of_route),\
        "w+")
    new_file.write("name of route: {}\n\n".format(name_of_route))
    new_file.write("origin coordinates: {}\n\n".format(origin_coordinates_input))
    new_file.write("destination coordinates: {}\n\n".format(\
        destination_coordinates_input))
    new_file.write("{}\n\n".format(dataset_misc_values))
    
    new_file.write("shortest_route\n")
    new_file.write("number: {}\n".format(index_real_shortest_route))
    new_file.write("time (s): {}\n".format(shortest_route_total_time))
    new_file.write("distance (m): {}\n".format(shortest_route_distance))
    new_file.write("MISC mean (-): {}\n".format(shortest_route_misc_mean))
    new_file.write("MISC max (-): {}\n\n".format(shortest_route_misc_max))
    
    new_file.write("least_sickening_route\n")
    new_file.write("number: {}\n".format(index_least_sickening_route))
    new_file.write("time (s): {}\n".format(least_sickening_route_total_time))
    new_file.write("distance (m): {}\n".format(least_sickening_route_distance))
    new_file.write("MISC mean (-): {}\n".format(least_sickening_route_misc_mean))
    new_file.write("MISC max (-): {}\n\n".format(least_sickening_route_misc_max))
    
    new_file.write("minimum_misc_max_route\n")
    new_file.write("number: {}\n".format(index_minimum_misc_max_route))
    new_file.write("time (s): {}\n".format(minimum_misc_max_route_total_time))
    new_file.write("distance (m): {}\n".format(minimum_misc_max_route_distance))
    new_file.write("MISC mean (-): {}\n".format(minimum_misc_max_route_misc_mean))
    new_file.write("MISC max (-): {}\n".format(minimum_misc_max_route_misc_max))
    new_file.close()
    
    dataset_end.to_pickle("{}/dataset_end.pkl".format(dirName))
    dataset_results.to_pickle("{}/dataset_results.pkl".format(dirName))
    
    return(dataset_results, dataset_end)