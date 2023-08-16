# -*- coding: utf-8 -*-
"""
Created on Mon Mar 20 13:40:01 2023

@author: oscar
"""

import os
import sys
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import flopy
from flopy.utils.triangle import Triangle as trng
import geopandas as gpd
from shapely.geometry import Polygon
import sfrmaker
import pandas as pd
import dsvFunctions as dsvf

def triangleMesh_to_Shapefile(cells2d,vertices,fileName):
    
    triangles= list()
    for cell in cells2d:
        vy = [vertices[cell[-1]][-1],vertices[cell[-2]][-1],vertices[cell[-3]][-1],vertices[cell[-1]][-1]]
        vx = [vertices[cell[-1]][-2],vertices[cell[-2]][-2],vertices[cell[-3]][-2],vertices[cell[-1]][-2]]
        polygon = Polygon(list(zip(vx,vy)))
        triangles.append(polygon)
    poly_df = gpd.GeoDataFrame(geometry=triangles)
    poly_df.to_file(fileName,driver='ESRI Shapefile')
    
def triangleMesh_to_Polygon(cells2d,vertices,fileName):
    triangles= list()
    for cell in cells2d:
        vy = [vertices[cell[-1]][-1],vertices[cell[-2]][-1],vertices[cell[-3]][-1],vertices[cell[-1]][-1]]
        vx = [vertices[cell[-1]][-2],vertices[cell[-2]][-2],vertices[cell[-3]][-2],vertices[cell[-1]][-2]]
        polygon = Polygon(list(zip(vx,vy)))
        triangles.append(polygon)
    return triangles
        
def create_Mesh_Dataframe(grid):
    grid.get_vertices()
    cell2d = grid.get_cell2d()
    vertices = grid.get_vertices()
    triangles= list()
    cell_ID = list()
    i=0
    for cell in cell2d:
        
        vy = [vertices[cell[-1]][-1],vertices[cell[-2]][-1],vertices[cell[-3]][-1],vertices[cell[-1]][-1]]
        vx = [vertices[cell[-1]][-2],vertices[cell[-2]][-2],vertices[cell[-3]][-2],vertices[cell[-1]][-2]]
        polygon = Polygon(list(zip(vx,vy)))
        triangles.append(polygon)
        cell_ID.append(i)
        i = i + 1
    grid_df = pd.DataFrame({"nodenumber":cell_ID,"geometry":triangles})
    return grid_df
    
def geoJson_to_ShapeFile(geoJsonPath):
    gdf = gpd.read_file(geoJsonPath)
    gdf.to_file('NonCordShapefile/file.shp')



geoJsonPath = "D:/Master_Erasmus/Thesis/example_basin/river_network.geojson"
mesh = "D:/Master_Erasmus/Thesis/example_basin/mesh.vtu"
geoJson_to_ShapeFile(geoJsonPath)
grid = dsvf.vtk_to_disv(mesh,cell_type='triangle',z_interpretation='top')
fileName= "NonCordShapefile/gridv2"
triangleMesh_to_Shapefile(grid["cells2d"],grid["vertices"],fileName)


