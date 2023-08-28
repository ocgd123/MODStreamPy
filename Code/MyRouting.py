# -*- coding: utf-8 -*-
"""
Created on Wed May 10 12:47:14 2023

@author: oscar
"""
import warnings
import os
import numpy as np
import flopy
import geopandas as gpd
from shapely.geometry import Polygon
import pandas as pd
import Routing_Classes as rt_Classes
import River_conections as rvrCon
import Routing_Calculation as rtCalc
import dsvFunctions as dsvf
from flopy.utils.triangle import Triangle as trng
import pyvista as pv
import shutil
import yaml

    
def create_GW_Backup(model_path,ws):
    """
    Under development 
    
    Create a backup of a groundwater model directory.
    
    This function creates a backup of a groundwater model directory by copying its contents to a new directory
    within the specified workspace.
    
    Parameters:
        model_path (str): Path to the groundwater model directory to be backed up.
        ws (str): Workspace where the backup directory will be created.
    
    Returns:
        None
    """
    ws_dir = [name for name in os.listdir(ws) if os.path.isdir(name)]
    if model_path in ws_dir:
        destination_dir = ""
        os.mkdir()
        shutil.copytree(model_path, destination_dir)
    
def get_boundary_cellid(gwf):
    """
    Retrieve the cell IDs of head boundary cells from a MODFLOW groundwater flow model.
    
    This function extracts the cell IDs of boundary cells from the specified MODFLOW groundwater flow model.
    
    Parameters:
        gwf (flopy.modflow.Modflow): MODFLOW groundwater flow model object.
    
    Returns:
        list: A list containing the cell IDs of boundary cells in the model.
    """
    boundsh = gwf.get_package("chd_0")
    bounds_cells = boundsh._data_list[-1].get_data()[0]
    cellids = []
    for cellid in bounds_cells:
        cellids.append(cellid[0][-1])
        print(cellid)
    return cellids


def run_Modflow_with_leak(river_Network,time_manager,stress_period, result_stor):
    """
    Run MODFLOW groundwater flow simulation with specified leak values.
    
    This function runs a MODFLOW groundwater flow simulation with the specified leak values from a river network.
    
    Parameters:
        river_Network (RiverNetwork): River network object containing information about river reaches.
        time_manager (TimeManager): Time manager object containing time-related information.
        stress_period (int): The stress period for which the simulation is run.
        result_stor (ResultsManager): Results manager object containing simulation results.
    
    Returns:
        None
    """
    rch_strs = {}
    qgwf_times = result_stor.results_dict[stress_period]["qgwf"]
    acum_flows = np.zeros(result_stor.num_cells)
    for qgwfs in time_manager.modflow_stperiods:
        acum_flows = acum_flows+qgwfs
        # for reach in river_Network.reaches.keys():
        #     cell = river_Network.reaches[reach].reach_cell
        #     qgwf = result_stor.qleak[0][cell]
        #     rch_data.append([0,cell,qgwf])
    
    #     rch_strs[stress_per] = rch_data
            
    # mfwell  = flopy.mf6.ModflowGwfwel(gwf,stress_period_data=rch_strs)
    # sim.write_simulation()
    # success, buff = sim.run_simulation()
    # gwf.remove_package('wel_0')
    # sim.write_simulation()
    
    rch_data = []
    for reach in river_Network.reaches.keys():
        qgwf = river_Network.reaches[reach].qgwf
        qgwf_converted = time_manager.convert_from_Seconds(time_manager.gw_units, qgwf)
        rch_data.append([0,river_Network.reaches[reach].reach_cell,qgwf_converted])
        
    rch_strs = {0:rch_data}
    mfwell  = flopy.mf6.ModflowGwfwel(gwf,stress_period_data=rch_strs,auto_flow_reduce=0)
    sim.write_simulation()
    success, buff = sim.run_simulation()
    gwf.remove_package('wel_0')
    sim.write_simulation()

def read_YAML(path):
    
    with open(path, "r") as stream:
        try:
            files = yaml.safe_load(stream)
            return files     
        
        except yaml.YAMLError as exc:
            print(exc)
'''
Creation of the classes based on the data frame and the connection dict
'''
crs = 'epsg:25833'
prjct_name = "Simple_Case_dummy"

yaml_files_path = r"D:/Master_Erasmus/Thesis/Simple_case/preprocess_files/simple_case_input_files.yaml"
files = read_YAML(yaml_files_path)

riv_hydr_prop_path = r"D:/Master_Erasmus/Thesis/Simple_case/preprocess_files/simple_case_river_inputs.yaml"
riv_hydr_prop = read_YAML(riv_hydr_prop_path)

#River boundaries Inflows
# inflow_dict = {0:50, 6:50}
inflow_dict = {0:10}

# File reading

river_Conections = files["river_sfrmaker_shapefile"]
sfr_df_path = files["river_sfrmaker_shapefile_csv"]
mesh_vtk_path = files["mesh_vtk_path"]
modflow_model_path = files["modflow_model_path"] 
exe_path  = files["exe_path"] 
result_path = files["result_path"]

# result_path = "D:/Master_Erasmus/Thesis/Simple_case/Results"
if files["diversion"] == "No_diversion_file_created":
    diversion_df = pd.DataFrame()
    
if files["diversion"] != "No_diversion_file_created":
    diversion_df = pd.read_csv(files["diversion"])



# #Synthetic Network
# crs = 'epsg:25833'
# prjct_name = "Synthetic_Network_Conductivity_width_0.1m_dummy"
# river_Conections = "D:/Master_Erasmus/Thesis/preprocess_files/mf_river.shp"
# sfr_df_path = "D:/Master_Erasmus/Thesis/preprocess_files/mf_df.csv"
# diversion_df_path = "D:/Master_Erasmus/Thesis/preprocess_files/mf_div_df.csv"
# mesh_vtk_path = r"D:/Master_Erasmus/Thesis/example_basin/mesh.vtu"
# modflow_model_path = r'D:/Master_Erasmus/Thesis/workSpaceRouting'
# exe_path  = "D:/mf6.4.1/bin/mf6.exe"
# result_path = "D:/Master_Erasmus/Thesis/Results"
# diversion_df = pd.read_csv(diversion_df_path)

'''
Ground water model
'''

#Call modflow model
#To do: Create a backup of the model in a temp file
sim = flopy.mf6.MFSimulation.load(sim_ws=modflow_model_path,exe_name=exe_path)
ic_list = []
gwf = sim.get_model()
sim.write_simulation()
success, buff = sim.run_simulation()
ic_list.append(sim.get_model().output.head().get_data()[0][0])
bounds_cells = get_boundary_cellid(gwf)

'''
Time Manager
'''
time_manager= rt_Classes.time_Manager(sim = sim)
stress = time_manager.modflow_stperiods
'''
Routing Creation
'''

sfr_df = pd.read_csv(sfr_df_path)
conect_data = rvrCon.river_Conections(river_Conections,crs = crs)

#Model characteristics
manning = [0.04 for i in range(0,len(sfr_df["name"]))]
geometry = ["rectangle" for i in range(0,len(sfr_df["name"]))]
conductivity_list = [0.0001 for i in range(0,len(sfr_df["name"]))]
rb_elev_list = [2 for i in range(0,len(sfr_df["name"]))]
thickness = [0.5 for i in range(0,len(sfr_df["name"]))]
slopes = [0.01 for i in range(0,len(sfr_df["name"]))]
width = [0.1 for i in range(0,len(sfr_df["name"]))]
conducivity_aquifer_list = [0.0001 for i in range(0,len(sfr_df["name"]))]
top_list = [10 for i in range(0,len(sfr_df["name"]))]

river_Network = rt_Classes.river_Network()

river_Network.create_Network(name_list= sfr_df["name"], slope_list = sfr_df["slope"], 
                            reach_id_list = sfr_df["rno"]-1, 
                            reach_cell_list = sfr_df["node"], 
                            length_list=sfr_df["rchlen"],
                            width_list = sfr_df["width"],
                            rbd_elev_list=sfr_df["strtop"],
                            conect_data =conect_data,riv_hydr_prop=riv_hydr_prop,
                            diversion_df = diversion_df)


'''
Result Storage

'''
mesh = pv.read(mesh_vtk_path)
num_cells = mesh.GetNumberOfCells()
result_stor = rt_Classes.results_Manager(time_manager.modflow_stperiods,
                                          time_manager.modflow_times,num_cells,result_path)

'''
River model flow inputs
'''


river_Network.add_Boundaries(inflow_dict)
river_Network.add_Reach_in_Boundaries(bounds_cells)
river_Network = rtCalc.calculate_Network(river_Network, conect_data)


hriv_i = np.zeros(mesh.GetNumberOfCells())
for reach in river_Network.reaches.keys():
    hriv_i[river_Network.reaches[reach].reach_cell] = river_Network.reaches[reach].hriv



'''
Write inflow from river
'''
hgwi = sim.get_model().output.head().get_data()[0][0]
    
'''
Coupled calculation
'''

# while iters < 15 and river_Network.rsme_qgwf>10**-6 and river_Network.rsme_hgw>10**-6 and river_Network.rsme_hriv>10**-6:
for stress_period in time_manager.modflow_stperiods:
    iters = 0
    stress_time = time_manager.modflow_times[stress_period][-1]
    while iters < 39 :
        time_id = 0
        for time in time_manager.modflow_times[stress_period]:
            
            river_Network.update_gw_Heads(sim,time)  
            result_stor.results_dict_prev[stress_period]["hgw"][time] = result_stor.results_dict[stress_period]["hgw"][time].copy()
            result_stor.results_dict[stress_period]["hgw"][time] = sim.get_model().output.head().get_data(totim = time)[0][0] 


            result_stor.results_dict_prev[stress_period]["hriv"][time] = result_stor.results_dict[stress_period]["hriv"][time]
            result_stor.results_dict_prev[stress_period]["qout"][time] = result_stor.results_dict[stress_period]["qout"][time]
            result_stor.results_dict_prev[stress_period]["qgwf"][time] = result_stor.results_dict[stress_period]["qgwf"][time]
            result_stor.results_dict_prev[stress_period]["qin"][time] = result_stor.results_dict[stress_period]["qin"][time]
            river_Network = rtCalc.calculate_Network(river_Network, conect_data)
            
            qgwfs = np.zeros(mesh.GetNumberOfCells())
            hriv = np.zeros(mesh.GetNumberOfCells())
            qout = np.zeros(mesh.GetNumberOfCells())
            qin = np.zeros(mesh.GetNumberOfCells())
            
            #Update result storage 
            for reach in river_Network.reaches.keys():
                
                qgwfs[river_Network.reaches[reach].reach_cell] = -(river_Network.reaches[reach].qgwf+qgwfs[river_Network.reaches[reach].reach_cell])
                hriv[river_Network.reaches[reach].reach_cell] = river_Network.reaches[reach].hriv + river_Network.reaches[reach].rbd_elev
                qout[river_Network.reaches[reach].reach_cell] = river_Network.reaches[reach].qout
                qin[river_Network.reaches[reach].reach_cell] = sum(river_Network.reaches[reach].qin)
                
            result_stor.results_dict[stress_period]["hriv"][time] = hriv
            result_stor.results_dict[stress_period]["qout"][time] = qout
            result_stor.results_dict[stress_period]["qgwf"][time] = qgwfs
            result_stor.results_dict[stress_period]["qin"][time] = qin

            time_id = time_id + 1
        
        
        run_Modflow_with_leak(river_Network,time_manager,stress_period,result_stor)
        river_Network.update_gw_Heads(sim,stress_time)    

    
        simulation_result = mesh.copy()

        
        iters = iters + 1
    print("Number of iters" + str(iters))
    # river_Network.restart_Network()
    
result_stor.export_Data_csv(prjct_name,river_Network)
result_stor.export_Data_vtu(mesh,prjct_name)

    
river_Network.calculate_Water_Balance()
# result_df = pd.DataFrame( {"gw_head (m)":head.tolist(),
#                           "h_stage (m)":hstage.tolist(),
#                           "q_leak_(m3/s)":qgwfs.tolist()})

# result_df.to_csv(os.path.join( result_path,'sim_results.'+prjct_name+'.csv'))
# for reach in river_Network.reaches.keys():
#     q_gwt = q_gwt + river_Network.reaches[reach].qgwf

# model = sim.get_model()
# grid = sim.get_package("disv")
# time = sim.get_package("tdis").perioddata.get_file_entry()
# units = sim.get_package("tdis")
# sim.simulation_data.mfdata[name, 'HDS', 'HEAD']

# '''
# MODFLOW model set up
# '''
# #Creation of the modflow test case
# name = 'mf'
# exe_path  = "D:/mf6.4.1/bin/mf6.exe"
# sim = flopy.mf6.MFSimulation(sim_name=name, version='mf6',
#                               exe_name=exe_path,
#                               sim_ws="D:/Master_Erasmus/Thesis/workSpaceRouting")

# tdis = flopy.mf6.ModflowTdis(sim, time_units='SECONDS',
#                               perioddata=[[1.0, 1, 1.]])

# gwf = flopy.mf6.ModflowGwf(sim, modelname=name, save_flows=True)

# ims = flopy.mf6.ModflowIms(sim, print_option='SUMMARY', complexity='complex', 
#                             outer_hclose=1.e-8, inner_hclose=1.e-8)

# tri = dsvf.vtk_to_disv(r"D:/Master_Erasmus/Thesis/example_basin/mesh.vtu")

# cell2d = tri["cells2d"] 

# #To do: make te clock whise polygon general
# for cell in range(0,len(cell2d)):
#     k  = cell2d[cell][-1]
#     l = cell2d[cell][-2]
#     cell2d[cell][-1] = l
#     cell2d[cell][-2] = k

# vertices = tri["vertices"]
# nvert = tri["nvert"]
# ncpl = tri["ncpl"]
# dis = flopy.mf6.ModflowGwfdisv(gwf, length_units='METERS',
#                                 nlay=1, ncpl=ncpl, nvert=nvert,
#                                 top=[5], botm=[0], 
#                                 vertices=vertices, cell2d=cell2d)

# npf = flopy.mf6.ModflowGwfnpf(gwf, k=0.0001)
# ic = flopy.mf6.ModflowGwfic(gwf, strt=10.0)
# sto = flopy.mf6.ModflowGwfsto(gwf,iconvert =1,ss = 0.)

# chdList = []

# up = [0,20,40,60,80]
# down = [19,39,59,79,99]


# for icpl in up: chdList.append([(0, icpl), 1])
# for icpl in down: chdList.append([(0, icpl), 0])

# chd = flopy.mf6.ModflowGwfchd(gwf, stress_period_data=chdList)
# oc = flopy.mf6.ModflowGwfoc(gwf,
#                             budget_filerecord='{}.cbc'.format(name),
#                             head_filerecord='{}.hds'.format(name),
#                             saverecord=[('HEAD', 'all'),
#                                         ('BUDGET', 'all')],
#                             printrecord=[('HEAD', 'all'),                  
#                                           ('BUDGET', 'all')])
# '''
# Simulaion set up
# '''
# sim.write_simulation()
# success, buff = sim.run_simulation()

