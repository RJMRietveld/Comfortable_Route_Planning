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
import numpy as np

"""
Run definitions of this module
"""
def run_route_generator_definitions(origin_coordinates_input,\
    destination_coordinates_input):
    
    coordinates_midpoint, distance_half, bearing = calculate_midpoint(\
        origin_coordinates_input, destination_coordinates_input)
        
    G, G_4, origin_node, destination_node, bbox_coordinates_nw,\
        bbox_coordinates_se = create_graph(origin_coordinates_input,\
        destination_coordinates_input, distance_half)
    
    routes =\
        create_routes(G, origin_coordinates_input, destination_coordinates_input,\
        origin_node, destination_node, coordinates_midpoint, distance_half, bearing)
    
    G = project(G)
    
    return(G, G_4, routes)

"""
Definitions of this module
"""
def calculate_midpoint(origin_coordinates_input, destination_coordinates_input):
    """
    Determines the coordinates of the midpoint, the half straight-line distance
    from the origin to the destination and the bearing/angle between the north-
    line and the straight-line from the origin to the destination.

    Parameters
    -------
    origin_coordinates_input : tuple
        Latitudal and longitudal coordinates of the origin.
    destination_coordinates_input : tuple
        Latitudal and longitudal coordinates of the destination.
    
    Returns
    -------
    coordinates_midpoint : tuple
        Latitudal and longitudal coordinates of the midpoint.
    distance_half : 
        The half straight-line distance from the origin to the destination in
        meters.
    bearing : 
        The angle between the north-line and the straight-line from the 
        origin to the destination in degrees.
    """
    
    radius_of_earth = 6378.1
    
    origin_lat = np.radians(origin_coordinates_input[0])
    origin_lon = np.radians(origin_coordinates_input[1])
    destination_lat = np.radians(destination_coordinates_input[0])
    destination_lon = np.radians(destination_coordinates_input[1])
    
    delta_lon = destination_lon - origin_lon
    delta_lat = destination_lat - origin_lat
    
    a = np.sin(delta_lat / 2)**2 + np.cos(origin_lat) * np.cos(destination_lat)\
        * np.sin(delta_lon / 2)**2
    c = 2 * np.arctan2(np.sqrt(a),np.sqrt(1 - a))
    
    distance_origin_to_destination = radius_of_earth * c
    
    distance_half = distance_origin_to_destination/2
    
    bearing = ox.bearing.get_bearing(origin_coordinates_input,\
        destination_coordinates_input)
    
    bearing_input = [bearing,bearing]
    
    coordinates_midpoint = calculate_coordinates_from_distance(\
        origin_coordinates_input, distance_half, bearing_input)
    
    return(coordinates_midpoint, distance_half, bearing)

def calculate_coordinates_from_distance(coordinates_input, distance, bearing):
    """
    Get coordinates of a location based on the coordinates of a reference point,
    the distance from this point and the bearing between the reference point and 
    the desired location.

    Parameters
    -------
    coordinates_input : tuple
        Coordinates of the reference point.
    distance : 
        The distance from the reference point to the new location in meters.
    bearing :
        The angle between the reference point and the desired location in degrees.
        
    Returns
    -------
    coordinates_output : tuple
        Coordinates of the new location.
    """
    
    radius_of_earth = 6378.1
    
    bearing_lat = np.pi*bearing[0]/180
    bearing_lon = np.pi*bearing[1]/180
    
    coordinates_input_lat = np.radians(coordinates_input[0]) 
    coordinates_input_lon = np.radians(coordinates_input[1])
    
    coordinates_output_lat = np.arcsin(np.sin(coordinates_input_lat)*np.cos(\
        distance/radius_of_earth) + np.cos(coordinates_input_lat)*np.sin(\
        distance/radius_of_earth)*np.cos(bearing_lat))
    
    coordinates_output_lon = coordinates_input_lon + np.arctan2(np.sin(\
        bearing_lon)*np.sin(distance/radius_of_earth)*np.cos(coordinates_input_lat),\
        np.cos(distance/radius_of_earth)-np.sin(coordinates_input_lat)*np.sin(\
        coordinates_output_lat))
    
    coordinates_output = np.degrees(coordinates_output_lat), np.degrees(\
        coordinates_output_lon)

    return(coordinates_output)
    
def create_bbox(origin_coordinates_input, destination_coordinates_input,\
    distance_half):
    """
    Get coordinates of the north-west and south-east corners of the Bounded
    Box (BBOX). The coordinates are determined on a certain distance from the
    origin and destination coordinates.

    Parameters
    -------
    origin_coordinates_input : tuple
        Latitudal and longitudal coordinates of the origin.
    destination_coordinates_input : tuple
        Latitudal and longitudal coordinates of the destination.
    distance_half : 
        Distance between the origin and destination coordinates and the BBOX.
    
    Returns
    -------
    bbox_coordinates_nw : tuple 
        Latitudal and longitudal coordinates of the north-west corner of the BBOX.
    bbox_coordinates_se: tuple
        Latitudal and longitudal coordinates of the south-east corner of the BBOX.
    """
    
    bbox_coordinate_n = max(origin_coordinates_input[0],\
        destination_coordinates_input[0])
    bbox_coordinate_e = max(origin_coordinates_input[1],\
        destination_coordinates_input[1])
    bbox_coordinate_s = min(origin_coordinates_input[0],\
        destination_coordinates_input[0])
    bbox_coordinate_w = min(origin_coordinates_input[1],\
        destination_coordinates_input[1])
    
    bbox_coordinates_nw_input = [bbox_coordinate_n,bbox_coordinate_w]
    bbox_coordinates_se_input = [bbox_coordinate_s,bbox_coordinate_e]
    
    bbox_distance = distance_half*0.6
    
    #NW
    bbox_bearing = [0,270]
    bbox_coordinates_nw = calculate_coordinates_from_distance(\
        bbox_coordinates_nw_input, bbox_distance, bbox_bearing)
  
    #SE
    bbox_bearing = [180,90]
    bbox_coordinates_se = calculate_coordinates_from_distance(\
        bbox_coordinates_se_input, bbox_distance, bbox_bearing)
    
    return(bbox_coordinates_nw, bbox_coordinates_se)

def create_graph(origin_coordinates_input, destination_coordinates_input,\
    distance_half):
    """
    Creates a graph of the BBOX-area and the route from the origin to destination.
    The data of all nodes and edges from within this area are restracted from
    OpenStreetMap (OSM).

    Parameters
    -------
    origin_coordinates_input : tuple
        Latitudal and longitudal coordinates of the origin.
    destination_coordinates_input : tuple
        Latitudal and longitudal coordinates of the destination.
    distance_half : 
        The half straight-line distance from the origin to the destination in
        meters.
        
    Returns
    -------
    G : networkx.MultiDiGraph
        Input graph.
    G_4 :  networkx.MultiDiGraph
        Input graph without any transformations.
    origin_node :  
        Latitudal and longitudal coordinates and node id of the origin node.
    destination_node : 
        Latitudal and longitudal coordinates and node id of the destination node.
    bbox_coordinates_nw : tuple 
        Latitudal and longitudal coordinates of the north-west corner of the BBOX.
    bbox_coordinates_se: tuple
        Latitudal and longitudal coordinates of the south-east corner of the BBOX.
    """
    
    bbox_coordinates_nw, bbox_coordinates_se = create_bbox(\
        origin_coordinates_input, destination_coordinates_input, distance_half)
    
    G_1 = ox.graph_from_bbox(bbox_coordinates_nw[0], bbox_coordinates_se[0],\
        bbox_coordinates_nw[1], bbox_coordinates_se[1], network_type='drive',\
        simplify=False)  
    G_2 = ox.add_edge_speeds(G_1)
    G_3 = ox.add_edge_bearings(G_2)
    G_4 = ox.add_edge_travel_times(G_3)
    G = G_4
      
    origin_node = ox.get_nearest_node(G, origin_coordinates_input)
    destination_node = ox.get_nearest_node(G, destination_coordinates_input)
    
    return(G, G_4, origin_node, destination_node, bbox_coordinates_nw,\
           bbox_coordinates_se)

def create_routes(G, origin_coordinates_input, destination_coordinates_input,\
    origin_node, destination_node, coordinates_midpoint, distance_half, bearing): 
    """
    Create a number of different routes from origin to destination coordinate.
    The differentiation is done by making use of intermediate points wich are 
    randomly chosen at a lineperpendicular  to  the  straight-line  at  the  
    half-distance of the straight-line from the origin to destination by making 
    use of a normal (Gaussian) distribution.

    Parameters
    -------
    G : networkx.MultiDiGraph
        Input graph.
    origin_coordinates_input : tuple
        Latitudal and longitudal coordinates of the origin.
    destination_coordinates_input : tuple
        Latitudal and longitudal coordinates of the destination.
    origin_node :  
        Latitudal and longitudal coordinates and node id of the origin node.
    destination_node : 
        Latitudal and longitudal coordinates and node id of the destination node.
    coordinates_midpoint : tuple
        Latitudal and longitudal coordinates of the midpoint.
    distance_half : 
        The half straight-line distance from the origin to the destination in
        meters.
    bearing :
        The angle between the reference point and the desired location in degrees.
    
    Returns
    -------
    routes : list 
        List of routes.
    """
    
    number_of_routes = 30

    route_shortest = ox.shortest_path(G, origin_node, destination_node,\
        weight="travel_time")
    
    routes = []
    routes.append(route_shortest)
    
    for i in range(number_of_routes):
        try:
            route, intermediate_coordinate = route_via_midpoint(G,\
                origin_coordinates_input, destination_coordinates_input,\
                origin_node, destination_node, coordinates_midpoint,\
                distance_half, bearing)
            route = unique_route_maker(route)
        except Exception:
            print('>>> route_via_midpoint() raises an error')
            pass
        else:
            routes.append(route)
    
    return(routes)

def filter_duplicates(route):
    """
    Function for the filtering of routes for duplicates and unnecessary loops.   
   
    Parameters
    ----------
    route : list
        List of nodes of an unfiltered route. 
        
    Returns
    -------
    new_route : list 
        List of route nodes filtered of duplicates and unnecessary loops.
    """
    switch = 0
    route_appendable = []
    new_route = []
    check_route = route 
    while switch == 0:
        dub  = []
        for i in range(len(route)):
            if route[i] == check_route[-1]:           
                new_route = new_route + route
                switch = 1
            if route[i] in route[:i]:              
                dub.append(route[i])       
                for j in range(len(route)):        
                    if route[i+j+1] not in route[:(i+j+1)]:
                        route_appendable = route_appendable + route[:(i+j+1)]
                        route = route[(i+j+1):]                   
                        break
                    else:
                        dub.append(route[i+j+1])
                route_part=route_appendable[:(route_appendable.index(dub[-1])+1)]
                route_appendable = []
                new_route = new_route + route_part               
                break             
    return(new_route)

def unique_route_maker(route):
    """
    Function for running the routes through a duplicate filter  
    
    Parameters
    ----------
    routes : list 
        List of unfiltered routes. 
        
    Returns
    -------
    filtered_routes : list
        List of filtered routes. 
    """
    
    route_set = set(route)
    max_runs = 4 
    for i in range(max_runs): 
        contains_duplicates = len(route) != len(route_set)
        if contains_duplicates == False: 
            filtered_routes = route
            break 
        else: 
            filtered_routes = filter_duplicates(route)
    return(filtered_routes)

def project(G):
    """
    Projection from a lat,long-coordinate system to an utm-coordinate system.
    
    Parameters
    -------
    G : networkx.MultiDiGraph
        Input graph (lat,long-coordinate system).
        
    Returns
    -------
    G : networkx.MultiDiGraph
        Input graph (utm-coordinate system).
    """
    
    G = ox.projection.project_graph(G)
    gdf_nodes,gdf_edges = ox.graph_to_gdfs(G)
    
    G = ox.graph_from_gdfs(gdf_nodes, gdf_edges, graph_attrs=G.graph)
    for u, v, data in G.edges(keys=False, data=True):
        assert 'geometry' in data
    
    return(G)

def route_via_midpoint(G, origin_coordinates_input, destination_coordinates_input,\
    origin_node, destination_node, coordinates_midpoint, distance_half, bearing):
    """
    the fastest route is created from the origin to the midpoint and from the  
    midpoint to the destination. These routes are appended to each other to 
    create a single route which passes through the intermediate point.  

    Parameters
    -------
    G : networkx.MultiDiGraph
        Input graph.
    origin_coordinates_input : tuple
        Latitudal and longitudal coordinates of the origin.
    destination_coordinates_input : tuple
        Latitudal and longitudal coordinates of the destination.
    origin_node :  
        Latitudal and longitudal coordinates and node id of the origin node.
    destination_node : 
        Latitudal and longitudal coordinates and node id of the destination node.
    coordinates_midpoint : tuple
        Latitudal and longitudal coordinates of the midpoint.
    distance_half : 
        The half straight-line distance from the origin to the destination in
        meters.
    bearing : 
        The angle between the north-line and the straight-line from the 
        origin to the destination in degrees.
    
    Returns
    -------
    route : list
        List of nodes of a route.
    intermediate_coordinate : 
        Latitudal and longitudal coordinates of the randomly generated
        intermediate point.
    """
    
    number_midpoints = 1
    origin_coordinates_midpoint = origin_coordinates_input
    origin_node_midpoint = origin_node
    route = []
    
    for i in range(number_midpoints+1):
        
        if i < (number_midpoints):
                        
            intermediate_bearing = ox.bearing.get_bearing(\
                origin_coordinates_midpoint, destination_coordinates_input)
            
            intermediate_bearing = [intermediate_bearing,intermediate_bearing]
            
            distance_new_point = ((2*distance_half)/(number_midpoints+1))/abs(\
                np.cos((intermediate_bearing[0]-bearing)*(np.pi/180)))
            
            coordinates_scatterline = calculate_coordinates_from_distance(\
                origin_coordinates_midpoint, distance_new_point,\
                intermediate_bearing)
            
            distance_from_midpoint = np.random.normal(loc=0.0, scale = 0.5) *\
                distance_half*0.5
            
            bearing_scatterline = bearing + 90
            bearing_scatterline = [bearing_scatterline,bearing_scatterline]
            
            intermediate_coordinate = calculate_coordinates_from_distance(\
                coordinates_scatterline, distance_from_midpoint, bearing_scatterline)
            
            intermediate_node = ox.get_nearest_node(G, intermediate_coordinate)
            
            route_origin_to_destination_midpoint = ox.shortest_path(G,\
                origin_node_midpoint, intermediate_node, weight="travel_time")
            
            del route_origin_to_destination_midpoint[-1] 
            
            route = [*route, *route_origin_to_destination_midpoint]
            
            origin_coordinates_midpoint = intermediate_coordinate
            origin_node_midpoint = intermediate_node
            
        else:
            route_origin_midpoint_to_destination = ox.shortest_path(G,\
                origin_node_midpoint, destination_node, weight="travel_time")

            route = [*route, *route_origin_midpoint_to_destination]
            
    return(route, intermediate_coordinate)

def plot_route(G, route):
    """
    Plot a route along a graph.
    
    Parameters
    -------
    G : networkx.MultiDiGraph
        Input graph.
    route : list
        Route as a list of node IDs.
        
    Returns
    -------
    figure, axis : tuple
        Matplotlib figure, axis.
    """
    
    figure, axis = ox.plot_graph_route(G, route)
    
    return(figure, axis)