# -*- coding: utf-8 -*-
"""
Created on Mon Apr  3 10:49:39 2023

@author: oscar
"""
import os
import geopandas as gpd

def river_Conections(shapefile_Path,crs):
    """
    Extract river connections from a shapefile.
    
    This function extracts river connections from a shapefile containing river geometries reaches and node coordinates.
    It generates a dictionary indicating which reaches are connected to each other upstream and downstream.
    
    Parameters:
        shapefile_Path (str): Path to the shapefile containing reach geometries and node coordinates.
        crs (str or dict): Coordinate reference system of the shapefile.
    
    Returns:
        dict: A dictionary representing river connections with 'in' and 'out' keys for each node.
            The 'in' key lists reach that flow into the specified reach.
            The 'out' key lists reach that the specified reach flows into.
    """
    conections_dict = {}
    shapefile = gpd.read_file(shapefile_Path,crs = crs)
    reaches = []
    
    for i in range(0,len(shapefile["rno"])):
        conections_dict[shapefile["rno"][i]] ={"in":[],"out":[]}
        
    for i in range(0,len(shapefile["rno"])):
        
        #To do: Vectorization using shapely
        node_id_i = shapefile["rno"][i]
        i_line_coord_upstrm = shapefile["geometry"][i].coords[0]
        i_line_coord_downstrm = shapefile["geometry"][i].coords[-1]
        num_bifur = 0

        for j in range(i+1,len(shapefile["rno"])):
            
            j_line_coord_upstrm = shapefile["geometry"][j].coords[0]
            j_line_coord_downstrm = shapefile["geometry"][j].coords[-1]
            distance_down = ((i_line_coord_downstrm[0]-j_line_coord_upstrm[0])**2+(i_line_coord_downstrm[1]-j_line_coord_upstrm[1])**2)**0.5
            distance_up = ((i_line_coord_upstrm[0]-j_line_coord_downstrm[0])**2+(i_line_coord_upstrm[1]-j_line_coord_downstrm[1])**2)**0.5
            node_id_j = shapefile["rno"][j]
            
            if distance_down<=1*10**-6:
                conections_dict[node_id_j]["in"].append(node_id_i-1)
                conections_dict[node_id_i]["out"].append(-(node_id_j-1))
                    
            if distance_up<=1*10**-6:
                conections_dict[node_id_j]["out"].append(node_id_i-1)
                conections_dict[node_id_i]["in"].append(-(node_id_j-1))

    return conections_dict



