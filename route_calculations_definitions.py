"""
This file is part of the Bachelor End Project (BEP)

Comfortable  Route  Planning

from mechanical engineering of the Delft university of technology.

Authors: Bas de Jager, Maurits Neele, Mike Rietveld and Julian waas
"""

"""
Import modules
"""
import osmnx as ox
import pandas as pd
import numpy as np
from scipy.ndimage import gaussian_filter1d
from scipy.interpolate import UnivariateSpline
from shapely.geometry import LineString

"""
Run definitions of this module
"""
def run_route_calculations_definitions(G, G_4, route):
    
    north_coordinates, east_coordinates, interpolation_step = first_interpolation(\
        G, route)
    
    north_coordinates, east_coordinates = gaussian_filter(north_coordinates,\
        east_coordinates)
    
    north_coordinates, east_coordinates, ss, S = second_interpolation(\
        north_coordinates, east_coordinates, interpolation_step)
    
    generated_interpolation_step, delta_generated_interpolation_step,\
        minimum_generated_interpolation_step, maximum_generated_interpolation_step\
        = check_interpolation_steps(north_coordinates, east_coordinates,\
        interpolation_step)
    
    K, radius = curvature_splines(north_coordinates, east_coordinates)
    
    max_speed_per_node = get_edge_data(G_4, route, north_coordinates,\
        east_coordinates, ss)
    
    dataset = create_dataset(north_coordinates, east_coordinates,\
        generated_interpolation_step, K, max_speed_per_node)
    
    return(generated_interpolation_step, K, max_speed_per_node, dataset, ss)

"""
Definitions of this module
"""
def first_interpolation(G, route):
    """
    Interpolation evenly spaced points along a LineString based on 
    geometry of the path.
    The spacing is approximate because the LineString's length may not be
    evenly divisible by it.

    Parameters
    -------
    G : networkx.MultiDiGraph
        Input graph.
    route : list
        Route as a list of node IDs.
    interpolation_step : float
        Interpolation step given in meters.
    
    Returns
    -------
    north_coordinates : array of float
        North_coordinates coordinates of new generated nodes.
    east_coordinates : array of float
        East_coordinates coordinates of new generated nodes.
    """
    
    interpolation_step = 0.5
    new_nodes_first_interpolation = np.empty(shape=(1,2))
    
    geometry = ox.utils_graph.get_route_edge_attributes(G, route, "geometry")
    
    for i in range(len(geometry)):   
        route_first_interpolation =  interpolate_points(geometry[i],\
            interpolation_step)
        route_first_interpolation = np.array((list(route_first_interpolation)))
        new_nodes_first_interpolation = np.append(new_nodes_first_interpolation,\
            route_first_interpolation, axis = 0)
    
    north_coordinates = new_nodes_first_interpolation[1:-1,0]
    east_coordinates = new_nodes_first_interpolation[1:-1,1]
    
    return(north_coordinates, east_coordinates, interpolation_step)

def interpolate_points(geometry, interpolation_step):
    """
    Interpolate evenly spaced points along a LineString.
    The spacing is approximate because the LineString's length may not be
    evenly divisible by it.
    
    Parameters
    ----------
    geom : shapely.geometry.LineString
        A LineString geometry.
    dist : float
        Spacing distance between interpolated points, in same units as `geom`.
        smaller values generate more points.
    
    Returns
    ------
    points : generator
        A generator of (x, y) tuples of the interpolated points' coordinates.
    """
    if isinstance(geometry, LineString):
        num_vert = max(round(geometry.length / interpolation_step), 1)
        for n in range(num_vert + 1):
            point = geometry.interpolate(n / num_vert, normalized=True)
            yield point.x, point.y
    else:
        raise TypeError(f"unhandled geometry type {geometry.geom_type}")
    
    return()

def gaussian_filter(north_coordinates, east_coordinates):
    """
    1-D Gaussian filter.
    
    Parameters
    -------
    sigma : scalar
        Standard deviation for Gaussian kernel.
    north_coordinates : array of float
        North_coordinates coordinates of nodes.
    east_coordinates : array of float
        East_coordinates coordinates of nodes.

    Returns
    -------
    north_coordinates : array of float
        North_coordinates coordinates of new generated nodes.
    east_coordinates : array of float
        East_coordinates coordinates of new generated nodes.
        
    Improvements
    -------
    Improvements can be made for the sigma input. This can be determined
    numerically for example. The minimum curvature radius is equal to 
    5 meters. The sigma can be calculated by taking the difference in the 
    smallest radius of the route and the minimum curvature radius of 5 
    meter. The difference can be the sigma factor. 
    """

    standard_deviation = 14
    
    north_coordinates = gaussian_filter1d(north_coordinates, standard_deviation)
    east_coordinates = gaussian_filter1d(east_coordinates, standard_deviation)
    
    return(north_coordinates, east_coordinates)

def second_interpolation(north_coordinates, east_coordinates, interpolation_step):
    """
    One-dimensional linear interpolation.
    Interpolation is done for the north_coordinates and east_coordinates
        coordinates seperately. 
    
    Parameters
    -------
    north_coordinates : array of float
        North_coordinates coordinates of nodes.
    east_coordinates : array of float
        East_coordinates coordinates of nodes.
    interpolation_step : float
        Interpolation step given in meters.

    Returns
    -------
    north_coordinates : array of float
        North_coordinates coordinates of new generated nodes.
    east_coordinates : array of float
        East_coordinates coordinates of new generated nodes.
    S : array of float
    ss : array of float
    """
    
    S = np.append(np.array(0),np.cumsum(np.sqrt(np.square(np.diff(\
        north_coordinates)) + np.square(np.diff(east_coordinates))))) 
    ss = np.arange(0,int(np.floor(S[-1])), interpolation_step)
    
    north_coordinates = np.interp(ss, S, north_coordinates)
    east_coordinates = np.interp(ss, S, east_coordinates)
   
    return(north_coordinates, east_coordinates, ss, S)

def check_interpolation_steps(north_coordinates, east_coordinates,\
    interpolation_step):
    """
    
    >>> Control function <<<
    
    Check the interpolation steps generated by the second_interpolation()
    When the generated interpolations differ too much in comparisation
    with the interpolation_step, the interpolation is not valid.
    
    Parameters
    -------
    north_coordinates : array of float
        North_coordinates coordinates of nodes.
    east_coordinates : array of float
        East_coordinates coordinates of nodes.
    interpolation_step : float
        Interpolation step given in meters.

    Returns
    -------
    ds : array of float
        The interpolation steps between two nodes after the second 
        interpolation.
    ds_delta : array of float
        The difference between interpolation steps between two nodes  
        after the second interpolation and the chosen interpolation_step.
    minimum_ds : float
        Smallest interpolation step generated by the second interpolation.
        This value needs to be close to the chosen interpolation_step. 
    maximum_ds : float
        Smallest interpolation step generated by the second interpolation.
        This value needs to be close to the chosen interpolation_step.
    """
    
    delta_north_coordinates = np.zeros(len(north_coordinates))
    delta_east_coordinates = np.zeros(len(delta_north_coordinates))
    generated_interpolation_step = np.zeros(len(delta_north_coordinates))
    delta_generated_interpolation_step = np.zeros(len(delta_north_coordinates))
    
    for i in range(len(delta_north_coordinates)-1):
        delta_north_coordinates[i] = north_coordinates[i+1]-north_coordinates[i]
        delta_east_coordinates[i] = east_coordinates[i+1]-east_coordinates[i]
        generated_interpolation_step[i] = np.sqrt((delta_north_coordinates[i])**2+\
            (delta_east_coordinates[i])**2)
        delta_generated_interpolation_step[i] = generated_interpolation_step[i]-\
            interpolation_step
        
        if round(generated_interpolation_step[i],2) != interpolation_step:
            print(">>> Interpolation_step is incorrect :",\
                  generated_interpolation_step[i])
            
    minimum_generated_interpolation_step = min(generated_interpolation_step[1:-1])
    maximum_generated_interpolation_step = max(generated_interpolation_step[1:-1])
    
    return(generated_interpolation_step, delta_generated_interpolation_step,\
           minimum_generated_interpolation_step,\
               maximum_generated_interpolation_step)

def curvature_splines(north_coordinates, east_coordinates):
    """
    Calculate the signed curvature of a 2D curve at each point
    using interpolating splines.
    
    Parameters
    ----------
    north_coordinates : array of float
        North_coordinates coordinates of nodes.
    east_coordinates : array of float
        East_coordinates coordinates of nodes.
    error : float
        The admisible error when interpolating the splines.
    w : (N,) array_like, optional
        Weights for spline fitting.  Must be positive.  If None (default),
        weights are all equal.
    k : int, optional
        Degree of the smoothing spline.  Must be 1 <= `k` <= 5.
        Default is `k` = 3, a cubic spline.
    
    Returns
    -------
    K : array of float
        Curvature (in meters^(-1)).
    radius : array of float
        Radius (in meters).
        
    Improvements
    -------
    The error is chosen based on certain tests, for example making use 
    of a circle as input for this function. The curvature values could 
    be validaded visually. An improvement can be made on the determination
    of the error value.
    """

    error=0.000001    
    
    t = np.arange(north_coordinates.shape[0])
    std = error * np.ones_like(north_coordinates)

    fx = UnivariateSpline(t, north_coordinates, k=4, w=1 / np.sqrt(std))
    fy = UnivariateSpline(t, east_coordinates, k=4, w=1 / np.sqrt(std))

    xˈ = fx.derivative(1)(t)
    xˈˈ = fx.derivative(2)(t)
    yˈ = fy.derivative(1)(t)
    yˈˈ = fy.derivative(2)(t)
    K = (xˈ* yˈˈ - yˈ* xˈˈ) / np.power(xˈ** 2 + yˈ** 2, 3 / 2)
    
    radius = 1/K
    
    minimum_radius = 5
    
    minimum_radius_route = min(abs(radius))
    
    if minimum_radius_route <= minimum_radius:
        print(">>> The minimum radius of the route =",minimum_radius_route)
    
    return(K, radius)

def get_edge_data(G_4, route, north_coordinates, east_coordinates, ss):
    """
    
    Warning
    ----------
    The function get_nearest_edges() takes a long time. This (less
    accurate) approach is used. 
    
    
    The edge speed is determined by making use the data of the edges 
    before interpolating and filtering. This is appended to the new nodes
    by making use of the distance from the start. 
    
    Parameters
    ----------
    G_4 : networkx.MultiDiGraph
        Input graph -> unmanipulated G.
    route : list
        Route as a list of node IDs.
    north_coordinates : array of float
        North_coordinates coordinates of nodes.
    S : array of float
    ss : array of float
    
    Returns
    -------
    max_speed_per_node : array of float
        Maximum speeds per node, while the graph_from_address() is 
        unsimplified, the speed_kph equals the max_speed. The data of
        speed_kph is more complete, so this data is used. 
    """
    
    edge_length = ox.utils_graph.get_route_edge_attributes(G_4, route, "length")
    edge_speed = ox.utils_graph.get_route_edge_attributes(G_4, route, "speed_kph")
    edge_length_cumsum = np.cumsum(edge_length)
    
    max_speed_per_node = np.zeros(len(north_coordinates)) 
    j = 0
    for i in range(len(north_coordinates)):
        max_speed_per_node[i] = edge_speed[j]
        if (ss[i])*edge_length_cumsum[-1]/ss[-1] >= edge_length_cumsum[j]:
            j = j+1
            
    max_speed_per_node = (max_speed_per_node)/(3.6)
    
    return(max_speed_per_node)

def create_dataset(north_coordinates, east_coordinates,\
    generated_interpolation_step, K, max_speed_per_node):
    """
    A dataset of the new edge information is created. This completes the
    route generation.
    
    Parameters
    ----------
    north_coordinates : array of float
        North_coordinates coordinates of nodes.
    east_coordinates : array of float
        East_coordinates coordinates of nodes.
    ds : array of float
        The interpolation steps between two nodes after the second 
        interpolation.
    K : array of float
        Curvature (in meters^(-1))
    max_speed_per_node : array of float
        Maximum speeds per node, while the graph_from_address() is 
        unsimplified, the speed_kph equals the max_speed. The data of
        speed_kph is more complete, so this data is used. 

    Returns
    -------
    dataset : Dataframe
        Dataframe which contains the values of the input parameters.
    """
    
    data = {"north_coordinates":north_coordinates,"east_coordinates":\
        east_coordinates,"generated_interpolation_step":\
        generated_interpolation_step,"K":K,"max_speed":max_speed_per_node}
    dataset = pd.DataFrame(data,columns=['north_coordinates','east_coordinates',\
        'generated_interpolation_step','K','max_speed'])
    
    return(dataset, north_coordinates, east_coordinates)