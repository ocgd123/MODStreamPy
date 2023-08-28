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
import yaml

def triangleMesh_to_Shapefile(tri,fileName,crs):
    """
    Convert a triangle mesh to a shapefile.
    
    This function converts a triangle mesh, represented by vertices and cell information,
    into a shapefile containing polygons that represent the triangles of the mesh.
    
    Parameters:
        tri (dict): A dictionary containing mesh information including vertices and cell2d.
        fileName (str): Path to the output shapefile.
        crs (str or dict): Coordinate reference system (CRS) of the shapefile.
    
    Returns:
        None
    """
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
    """
    Generate an empty MODFLOW model.
    
    This function creates an empty MODFLOW model with a single layer and a specified simulation workspace.
    
    Parameters:
        exe_path (str): Path to the MODFLOW executable.
        sim_ws (str): Path to the simulation workspace where model files will be saved.
    
    Returns:
        gwf (flopy.modflow.Modflow): An empty MODFLOW groundwater flow model object.
    """
    
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

def create_Connection_df(connection_data,div_csv_path,files_yaml):
    
    """
    Create a DataFrame describing connections between main and diversion reaches.
    
    This function takes connection data in the form of a dictionary and generates a DataFrame
    to describe the connections between main and diversion reaches. It identifies diversion
    reaches, their corresponding out_reaches, reach types, and diversion flow rates.
    
    Parameters:
        connection_data (dict): Dictionary containing connection information for each reach.
        div_csv_path (str): Path to save the generated CSV file.
    
    Returns:
        None
    """
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
        files_yaml[0]["diversion"] = div_csv_path
        
    return files_yaml
        
def geoJson_to_Pandas(geoJsonPath,geoJsonName):
    
    """
    Convert a GeoJSON file to a shapefile using GeoPandas.
    
    This function reads a GeoJSON file, converts it to a GeoDataFrame using GeoPandas, and then
    exports the GeoDataFrame as a shapefile.
    
    Parameters:
        geoJsonPath (str): Path to the directory containing the GeoJSON file.
        geoJsonName (str): Name of the GeoJSON file (without extension).
    
    Returns:
        None
    """
    gdf = gpd.read_file(geoJsonPath+chr(47)+geoJsonName)
    gdf.to_file(geoJsonPath+chr(47)+geoJsonName+'.shp')


'''
Data preprocessing
'''
def vtu_to_shapefile(ws,vtu_mesh_path,shp_name,crs):
    """
    Convert a VTK unstructured grid mesh to a shapefile.
    
    This function takes a VTK unstructured grid mesh in VTU format and converts it
    into a shapefile containing polygons that represent the mesh elements.
    
    Parameters:
        ws (str): Working directory where the shapefile will be saved.
        vtu_mesh_path (str): Path to the VTU unstructured grid mesh file.
        shp_name (str): Name of the shapefile to be created.
        crs (str): Coordinate Reference System (CRS) of the shapefile.
    
    Returns:
        None
    """
        
    os.listdir(ws)
    
    if  "preprocess_files" not in os.listdir(ws):
        file_path  = os.path.join(ws, "preprocess_files")
        os.mkdir(file_path)
    else:
        file_path = os.path.join(ws, "preprocess_files")   
        
    shp_path = os.path.join(file_path, shp_name)
    tri = dsvf.vtk_to_disv(vtu_mesh_path)
    triangleMesh_to_Shapefile(tri,shp_path,crs)

def create_YAML_input_files_dict():
    yaml_file = [
        {        
        "river_sfrmaker_shapefile":"Input",
        "mesh_vtk_path": "Input",
        "diversion" : "No_diversion_file_created",
        "river_sfrmaker_shapefile_csv":"Input",
        "modflow_model_path" : "Input",
        "exe_path": "Input",
        "result_path":"Input"
        }
        ]
    return yaml_file
    
def create_YAML_Reaches(river_df,preprocess_path,project_name):
    data = []
    for line in range(0,len(river_df)):
        print(line)
        reach_dict = {
        river_df["name"][line]:
            {
        "River properties":
            {
            "river_conductivity" : "input",
            "river_manning" : "input",
            "river_bed_thickness" : "input"
            }
            }
        }
        data.append(reach_dict)
    with open(os.path.join(preprocess_path, project_name+"_river_inputs.yaml"), 'w',) as reach_yaml :
        yaml.dump_all(data,reach_yaml, sort_keys=False)
#Input Modflow exe path
#Project directory (different from modflow ws)
def reach_assignation(exe_path,ws, grid_shapefile, river_shape_file, project_name, crs, vtu_mesh_path,write_sfr = False,):
    """
    Assign reaches to the grid cells and export the results.
    
    This function assigns river data from the provided shapefile to the grid cells and exports the
    assigned rivers as shapefiles and dataframes for further processing.
    
    Parameters:
        exe_path (str): Path to the MODFLOW executable.
        ws (str): Working directory where the files will be saved.
        grid_shapefile (str): Path to the shapefile containing grid cell information.
        river_shape_file (str): Path to the shapefile containing river data.
        project_name (str): Name of the project or dataset.
        crs (str): Coordinate Reference System (CRS) of the shapefiles.
        write_sfr (bool, optional): Flag to write SFR package data for MODFLOW. Default is False.
    
    Returns:
        dict: A dictionary containing reaches connections information.
    """
    
    if  "preprocess_files" not in os.listdir(ws):
        file_path = os.path.join(ws, "preprocess_files")
        os.mkdir(file_path)
    else:
        file_path = os.path.join(ws, "preprocess_files")
        
    '''
    Data preporcessing
    '''
    #YAML dictionary with the created files
    files_yaml = create_YAML_input_files_dict()
    
    #Shape file load
    grd_df2 = gpd.read_file(grid_shapefile)
    
    #River shapefile data_frame
    river_df = gpd.read_file(river_shape_file)
    
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
    
    files_yaml[0]["exe_path"] = exe_path
    files_yaml[0]["mesh_vtk_path"] = vtu_mesh_path
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
    
    #Export the assigned rivers to a shapefile file for post processing
    sfrdata = custom_lines.to_sfr(grid=grd, model=gwf)
    if write_sfr == True:
        sfrdata.write_package(version='mf6',filename='Test2')
        
    #Grid shapefile with only the cells that contain a reach
    post_process_grid = os.path.join(file_path, project_name+"_grid.shp")
    sfrdata.export_cells(filename=post_process_grid)
    
    #River shapefile with the river discretization into reaches
    post_process_river = os.path.join(file_path, project_name+"_river.shp")
    sfrdata.export_lines(filename = post_process_river)
    files_yaml[0]["river_sfrmaker_shapefile"] = os.path.join(file_path, project_name+"_river.shp")

    # Dataframe with the attribute table from the post_process_river shapefile
    post_process_df = os.path.join(file_path, project_name+"_df.csv")
    sfr_df = sfrdata.reach_data    
    sfr_df.to_csv(post_process_df)
    files_yaml[0]["river_sfrmaker_shapefile_csv"] = os.path.join(file_path, project_name+"_df.csv")
    conect_data = rvrCon.river_Conections(post_process_river, crs)
    
    div_csv_path = os.path.join(file_path, project_name+"_div_df.csv")
    files_yaml = create_Connection_df(conect_data,div_csv_path,files_yaml)
    
    with open(os.path.join(file_path, project_name+"_input_files.yaml"), 'w',) as files_inputs :
        yaml.dump_all(files_yaml,files_inputs, sort_keys=False)
        
    create_YAML_Reaches(river_df,file_path,project_name)
    return sfr_df

'''
Example
'''
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

connect = reach_assignation(modfl_exe_path,ws, grid_shap_path,river_shape_file, project_name,crs,vtu_mesh_path)

