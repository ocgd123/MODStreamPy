# -*- coding: utf-8 -*-
"""
Created on Sat May 20 09:28:27 2023

@author: oscar
"""
import os
import pandas as pd
import flopy
import geopandas as gpd
from shapely.geometry import Polygon
import sfrmaker
import River_conections as rvrCon
import dsvFunctions as dsvf

def triangleMesh_to_Shapefile(tri,fileName,crs):
    
    vertices = tri["vertices"]
    cell2d = tri["cells2d"] 
    vertices = tri["vertices"]
    triangles= list()

    for cell in cell2d:
        
        vy = [vertices[cell[-1]][-1],vertices[cell[-2]][-1],vertices[cell[-3]][-1],vertices[cell[-1]][-1]]
        vx = [vertices[cell[-1]][-2],vertices[cell[-2]][-2],vertices[cell[-3]][-2],vertices[cell[-1]][-2]]
        polygon = Polygon(list(zip(vx,vy)))
        triangles.append(polygon)
        
    poly_df = gpd.GeoDataFrame(geometry = triangles, crs = crs)
    poly_df.to_file(fileName,driver='ESRI Shapefile')

def generate_Empty_Modflow_Model(exe_path,sim_ws):
    
    name = 'mf'
    sim = flopy.mf6.MFSimulation(sim_name=name, version='mfusg',
                                 exe_name=exe_path,
                                 sim_ws=sim_ws)
    tdis = flopy.mf6.ModflowTdis(sim, time_units='SECONDS',
                                 perioddata=[[1.0, 1, 1.]])
    gwf = flopy.modflow.Modflow(modelname=name)
    gwf._version = 'mfusg'
    nlay = 1
    top = 1.
    botm = [0.]
    dis = flopy.modflow.ModflowDis(gwf,nlay=nlay,
                               top=top, botm=botm)
    bas = flopy.modflow.ModflowBas(gwf, ibound=1)
    return gwf

def create_Connection_df(connection_data,div_csv_path):
    div_reaches = []
    out_reaches = []
    #Reach type must be Main or diversion
    reach_type = []
    #div_q is 0 for main and can have any value for the diversion
    div_q = []
    for reach in connection_data.keys():
        if len(connection_data[reach]["out"])>1:
            for out_reach in connection_data[reach]["out"]:
                div_reaches.append(reach)
                out_reaches.append(abs(out_reach)+1)
                reach_type.append("main/diversion")
                div_q.append(0)
    
    if len(div_reaches)>0:
        df = pd.DataFrame()
        df["diverted_reach"] = div_reaches
        df["out_reach"] = out_reaches
        df["reach_type"] = reach_type
        df["q"] = div_q
        df.to_csv(div_csv_path)
    
def geoJson_to_Pandas(geoJsonPath,geoJsonName):
    gdf = gpd.read_file(geoJsonPath+chr(47)+geoJsonName)
    gdf.to_file(geoJsonPath+chr(47)+geoJsonName+'.shp')


'''
Data preprocessing
'''
def vtu_to_shapefile(ws,vtu_mesh_path,shp_name,crs):
    
    os.listdir(ws)
    
    if  "preprocess_files" not in os.listdir(ws):
        file_path  = os.path.join(ws, "preprocess_files")
        os.mkdir(file_path)
    else:
        file_path = os.path.join(ws, "preprocess_files")   
        
    shp_path = os.path.join(file_path, shp_name)
    tri = dsvf.vtk_to_disv(vtu_mesh_path)
    triangleMesh_to_Shapefile(tri,shp_path,crs)

#Input Modflow exe path
#Project directory (different from modflow ws)
def river_assignation(exe_path,ws, grid_shapefile, river_shape_file, project_name, crs, write_sfr = False):
    
    if  "preprocess_files" not in os.listdir(ws):
        file_path = os.path.join(ws, "preprocess_files")
        os.mkdir(file_path)
    else:
        file_path = os.path.join(ws, "preprocess_files")
        
    '''
    Data preporcessing
    '''
    
    #Shape file load
    grd_df2 = gpd.read_file(grid_shapefile)
    
    #Crete the MODFLOW model
    gwf = generate_Empty_Modflow_Model(modfl_exe_path,ws)
    
    #Cretes the unstructured grid from a shapefile to DISV
    grd = sfrmaker.grid.UnstructuredGrid.from_dataframe(
        grd_df2,
        node_col='nodenumber',
        geometry_column='geometry',
        model_units='meters',
        crs='epsg:25833'
    )
    
    #US shapeflie path  D:/Master_Erasmus/Thesis/sfrmaker/examples/meras/flowlines.shp
    #Synthetic river path D:/Master_Erasmus/Thesis/Germany_River1_v2.shp
    #Assignation of the rivers to the grid cells
    custom_lines = sfrmaker.Lines.from_shapefile(shapefile = river_shape_file,
                                                  id_column='COMID',  # arguments to sfrmaker.Lines.from_shapefile
                                                  routing_column='tocomid',
                                                  width1_column='width1',
                                                  width2_column='width2',
                                                  up_elevation_column='elevupsmo',
                                                  dn_elevation_column='elevdnsmo',
                                                  name_column='GNIS_NAME',
                                                  attr_length_units='meters',  # units of source data
                                                  attr_height_units='meters',
                                                  crs='epsg:25833',
                                                  )
    
    #m = flopy.modflow.Modflow.load('shellmound.nam', model_ws='D:/Master_Erasmus/Thesis/sfrmaker/sfrmaker/test/data/shellmound/shellmound')
    #Export the assigned rivers to a shapefile file for post processing
    sfrdata = custom_lines.to_sfr(grid=grd, model=gwf)
    if write_sfr == True:
        sfrdata.write_package(version='mf6',filename='Test2')
        
    post_process_grid = os.path.join(file_path, project_name+"_grid.shp")
    sfrdata.export_cells(filename=post_process_grid)
    
    post_process_river = os.path.join(file_path, project_name+"_river.shp")
    sfrdata.export_lines(filename = post_process_river)
    
    post_process_df = os.path.join(file_path, project_name+"_df.csv")
    sfr_df = sfrdata.reach_data    
    sfr_df.to_csv(post_process_df)
    
    conect_data = rvrCon.river_Conections(post_process_river, crs)
    
    div_csv_path = os.path.join(file_path, project_name+"_div_df.csv")
    create_Connection_df(conect_data,div_csv_path)
    return conect_data

ws ="D:/Master_Erasmus/Thesis/Simple_case"
crs = 25833
shp_name = "model_grid.shp"
vtu_mesh_path = r'D:/Master_Erasmus/Thesis/Simple_case/Simple_Grid.vtk'
vtu_to_shapefile(ws,vtu_mesh_path,shp_name,crs)

#Modflow EXE path
modfl_exe_path = "D:/mf6.4.1/bin/mf6.exe"
#Grid shape path
grid_shap_path = r'D:/Master_Erasmus/Thesis/Simple_case/Simple_Grid.shp'
#River shape file path
river_shape_file = r'D:/Master_Erasmus/Thesis/Simple_case/Simple_net.shp'
#sfrmaker requieres to include the espg in the name of the coordinate system
crs = 'epsg:25833'
project_name = "simple_case"

connect = river_assignation(modfl_exe_path,ws, grid_shap_path,river_shape_file, project_name,crs)

