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
import matplotlib.tri as mtri
import River_conections as rvrCon
import pyvista as pv




def triangleMesh_to_Shapefile(grid,fileName):
    
    grid.get_vertices()
    
    cell2d = grid.get_cell2d()
    
    vertices = grid.get_vertices()
    
    triangles= list()

    for cell in cell2d:
        
        vy = [vertices[cell[-1]][-1],vertices[cell[-2]][-1],vertices[cell[-3]][-1],vertices[cell[-1]][-1]]
        
        vx = [vertices[cell[-1]][-2],vertices[cell[-2]][-2],vertices[cell[-3]][-2],vertices[cell[-1]][-2]]
        
        polygon = Polygon(list(zip(vx,vy)))
        
        triangles.append(polygon)
        
    poly_df = gpd.GeoDataFrame(geometry=triangles)
    
    poly_df.to_file(fileName,driver='ESRI Shapefile')
        
def create_Mesh_Dataframe(grid):
    
    grid.get_vertices() 
    cell2d = grid.get_cell2d()
    vertices = grid.get_vertices()
    triangles= list()
    cell_ID = list()
    i  =list()
    k = list()
    nodei = 0
    nodek = 1
    
    for cell in cell2d:
        vy = [vertices[cell[-1]][-1],vertices[cell[-2]][-1],vertices[cell[-3]][-1],vertices[cell[-1]][-1]]
        vx = [vertices[cell[-1]][-2],vertices[cell[-2]][-2],vertices[cell[-3]][-2],vertices[cell[-1]][-2]]
        polygon = Polygon(list(zip(vx,vy)))
        triangles.append(polygon)
        cell_ID.append(nodei)
        k.append(nodei)
        i.append(nodei)
        nodei = nodei + 1
    
    grid_df = pd.DataFrame({"k":k,"i":i,"j":i,"nodenumber":cell_ID,"geometry":triangles})
    return grid_df
    
def geoJson_to_Pandas(geoJsonPath):
    
    gdf = gpd.read_file('file.geojson')
    
    gdf.to_file('file.shp')

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

def sfr_to_vtk(mesh,sfr,reach_cell_dict,project_name):
    
    simulation_result = mesh.copy()
    n_cells = mesh.n_cells
    record_names = sfr.output.budget().get_unique_record_names(decode = True)
    sfr_budget = sfr.output.budget()
    sfr_times = sfr_budget.get_times()
    river_cell_id = list(reach_cell_dict.values())
    count = 0
    bc_cells  = dict
    bc_cells.update({'river':{'bc_id':1,'cell_ids':river_cell_id}})
    init_time = True
    for time in sfr_times:
        if init_time == True:

            all_outputs  = sfr_budget.get_data(totim = time)
            sfr_heads = sfr.output.stage().get_data(totim = time)[0][0]

            for output in all_outputs:
                output_name = record_names[count]
                if "FLOW-JA-FACE" in output_name:
                    
                    outFlow  = np.zeros(n_cells)
                    inFlow  = np.zeros(n_cells)
                    flow = output["q"]
                    node2 = output["node2"]
                    for i in range(0,len(flow)):
                        if flow[i] <0:
                            outFlow[reach_cell_dict[node2[i]]] = flow[i] + outFlow[reach_cell_dict[node2[i]]]
                        if flow[i] >0:
                            inFlow[reach_cell_dict[node2[i]]] = flow[i] + inFlow[reach_cell_dict[node2[i]]] 
                    
                    simulation_result = mesh.copy()
                    simulation_result.cell_data[output_name]=outFlow
                    simulation_result.save(os.path.join("D:/Master_Erasmus/Thesis/Results",'sim_t_out_1.vtu'))
                    simulation_result = mesh.copy()
                    simulation_result.cell_data[output_name]=inFlow
                    simulation_result.save(os.path.join("D:/Master_Erasmus/Thesis/Results",'sim_t_in_1.vtu'))                    
                    
                if "GWF" in output_name:
                    output_array_gwf  = np.zeros(n_cells)
                    output_array_head  = np.zeros(n_cells)
                    bc_cells.update({output_name.lower():{'bc_id':count+1,'cell_ids':river_cell_id}})
                    flow = output["q"]
                    node = output["node"]
                    node2 = output["node2"]
                    for i in range(0,len(flow)):
                        output_array_gwf[node2[i]-1]=flow[i]
                        output_array_head[node2[i]-1]=sfr_heads[i]
                    simulation_result = mesh.copy()
                    simulation_result.cell_data["q_leak_SFR_(m3/s)"]=output_array_gwf
                    simulation_result.save(os.path.join("D:/Master_Erasmus/Thesis/Results",'sim_t_1_gwf.vtu'))
                    simulation_result = mesh.copy()
                    simulation_result.cell_data["h_stage_SFR (m)"]=output_array_head
                    simulation_result.save(os.path.join("D:/Master_Erasmus/Thesis/Results",'sim_t_1_head.vtu'))                    
                    print(output_array_gwf.tolist())
         
                else:
                    x=2
                    
                count = count + 1
                
            init_time = False

    # write_pvd(pvd_filename=os.path.join("D:/Master_Erasmus/Thesis/Results",'sim1.pvd'),ts_name_dict='sim_t_1.vtu')
    result_df = pd.DataFrame( {"gw_head (m)":sim.get_model().output.head().get_data()[0][0],
                             "h_stage (m)":output_array_head,
                             "q_leak_(m3/s)":output_array_gwf})
    result_df.to_csv(os.path.join("D:/Master_Erasmus/Thesis/Results",'MODFLOW_resultq0.1.csv'))

    sfr_package_convergence = sfr.output.package_convergence()
    sfr_obs = sfr.output.obs()

workspace ="D:/Master_Erasmus/Thesis/Meras"
'''
#US case

#Define the active domain
active_domain = [(330177, 1246872), (333018, 1246872), (333018, 1249483), (330177, 1249483)]

#Creates triangular grid
tri = trng(angle=20, model_ws=workspace, exe_name="D:/Master_Erasmus/Thesis/Meras/Exe/triangle.exe")
tri.add_polygon(active_domain)
tri.add_region((330277,1246972),0,maximum_area=1000) #coarse discretization
tri.build()

#Plot the grid
fig = plt.figure(figsize=(15, 15))
ax = plt.subplot(1, 1, 1, aspect='equal')
tri.plot(edgecolor='gray')
cell2d = tri.get_cell2d()
vertices = tri.get_vertices()

cell2d = tri.get_cell2d()
vertices = tri.get_vertices()
xcyc = tri.get_xcyc()
ncpl = tri.ncpl
nvert = tri.nvert

grd_df = create_Mesh_Dataframe(tri)
leftcells = tri.get_edge_cells(4) 
rightcells = tri.get_edge_cells(2)
'''
#Input exe path
exe_path = "D:/mf6.4.1/bin/mf6.exe"

#Input: Simulation workspace
sim_ws = "D:/Master_Erasmus/Thesis/ModflowModel"

#Input Coordiante system
crs = 'epsg:25833'

#Crete the MODFLOW model
gwf = generate_Empty_Modflow_Model(exe_path,sim_ws)



#triangleMesh_to_Shapefile(tri,"grid2.shp")
grd_df2 = gpd.read_file("D:/Master_Erasmus/Thesis/gridv2/gridv2.shp")

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
custom_lines = sfrmaker.Lines.from_shapefile(shapefile='D:/Master_Erasmus/Thesis/Germany_River1_v2.shp',
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
sfrdata.write_package(version='mf6',filename='Test2')
sfrdata.export_cells(filename="Test2_grid.shp")
sfrdata.export_lines(filename="Test2.shp")

sfr_df = sfrdata.reach_data


'''
MODFLOW model set up
'''
#Creation of the modflow test case
name = 'mf'
times = [1]
sim = flopy.mf6.MFSimulation(sim_name=name, version='mf6',
                             exe_name=exe_path,
                             sim_ws="D:/Master_Erasmus/Thesis/workSpace")

tdis = flopy.mf6.ModflowTdis(sim, time_units='SECONDS',
                             perioddata=[[1, 1, 1.]])

gwf = flopy.mf6.ModflowGwf(sim, modelname=name, save_flows=True)

ims = flopy.mf6.ModflowIms(sim, print_option='SUMMARY', complexity='complex', 
                           outer_hclose=1.e-8, inner_hclose=1.e-8)

tri = dsvf.vtk_to_disv(r"D:/Master_Erasmus/Thesis/example_basin/mesh.vtu")

cell2d = tri["cells2d"] 

#To do: make te clock whise polygon general
for cell in range(0,len(cell2d)):
    k  = cell2d[cell][-1]
    l = cell2d[cell][-2]
    cell2d[cell][-1] = l
    cell2d[cell][-2] = k

vertices = tri["vertices"]
nvert = tri["nvert"]
ncpl = tri["ncpl"]
dis = flopy.mf6.ModflowGwfdisv(gwf, length_units='METERS',
                                nlay=1, ncpl=ncpl, nvert=nvert,
                                top=[5], botm=[0], 
                                vertices=vertices, cell2d=cell2d)

npf = flopy.mf6.ModflowGwfnpf(gwf, k=0.0001)
ic = flopy.mf6.ModflowGwfic(gwf, strt=5.0)
sto = flopy.mf6.ModflowGwfsto(gwf,iconvert =1,ss = 0)

chdList = []

up = [0,20,40,60,80]
down = [19,39,59,79,99]


for icpl in up: chdList.append([(0, icpl), 1])
for icpl in down: chdList.append([(0, icpl), 0])

chd = flopy.mf6.ModflowGwfchd(gwf, stress_period_data=chdList)
oc = flopy.mf6.ModflowGwfoc(gwf,
                            budget_filerecord='{}.cbc'.format(name),
                            head_filerecord='{}.hds'.format(name),
                            saverecord=[('HEAD', 'all'),
                                        ('BUDGET', 'all')],
                            printrecord=[('HEAD', 'all'),                  
                                         ('BUDGET', 'all')])

'''
SFR set up
'''
connectionData = rvrCon.river_Conections("Test2.shp",crs = crs)
pckgData = []
pckgConections = []
pckgDiversion = []
reach_cell_dict = {}
num_cells = 100
reach_ID_list = ["No_reach"] *num_cells

for i in range(0,sfr_df.shape[0]):
    reach_ID_list[sfr_df["node"][i]] = sfr_df["rno"][i]

    rno = sfr_df["rno"][i]
    num_Connections = len(connectionData[rno]["in"])+len(connectionData[rno]["out"])
    
    data = (sfr_df["rno"][i]-1, 0, sfr_df["node"][i], sfr_df["rchlen"][i],
            0.1, 0.01,2,
            0.5,0.0001,0.04,num_Connections,1,0)
    reach_cell_dict[sfr_df["rno"][i]] = sfr_df["node"][i]
    
    if len(connectionData[rno]["out"]) > 1:
        
        data = (sfr_df["rno"][i]-1, 0, sfr_df["node"][i], sfr_df["rchlen"][i],
                0.1, 0.01,2,
                0.5,0.0001,0.04,num_Connections,1,1)
        
        pckgDiversion.append([sfr_df["rno"][i]-1, 0, 15,"UPTO"])
            
    if rno == 16:
        data = (sfr_df["rno"][i]-1, 0,sfr_df["node"][i], sfr_df["rchlen"][i],
                0.1, 0.01,2,
                0.5,0.0001,0.04,num_Connections,0,0)
    
         
    pckgData.append(data)
    pckgConections.append([rno-1,*connectionData[rno]["out"],*connectionData[rno]["in"]])
    
channel_flow_data = [[0,"inflow",50],[6,"inflow",50],[11,"diversion",0,50]]

sfr = flopy.mf6.ModflowGwfsfr(gwf,nreaches=sfr_df.shape[0],packagedata=pckgData,
                              connectiondata=pckgConections,diversions=pckgDiversion,
                              print_flows=True,print_input =True,print_stage=True,
                              save_flows=True,stage_filerecord="sfr_stage.hds",
                              budgetcsv_filerecord ='{}.csv'.format(name),
                              budget_filerecord =  '{}.cbc'.format(name+"_SFR"),
                              perioddata=channel_flow_data,package_convergence_filerecord='conv{}.csv'.format(name)
                              )

'''
Simulaion set up
'''
sim.write_simulation()
success, buff = sim.run_simulation()

hdobj = flopy.utils.HeadFile("D:/Master_Erasmus/Thesis/workSpace/mf.hds", precision='double')




'''
Simulaion set up
'''


sim.write_simulation()
success, buff = sim.run_simulation()


'''
Write Results
'''
# sfr_flows = flopy.utils.CellBudgetFile("D:/Master_Erasmus/Thesis/Simple_case/SimpleModflow/mf_SFR.cbc")
# flows = flopy.utils.CellBudgetFile("D:/Master_Erasmus/Thesis/Simple_case/SimpleModflow/mf_SFR.cbc").get_data(idx = 1)
# times = [10]
cell_id_list = range(1, 101)

result_dict = {}
prjct_name = "Synthetic_Network_Base_Case_SFR_width_0.1m"
stress_len = 1
result_path = "D:/Master_Erasmus/Thesis/Results_SFR"


csv_name = prjct_name+"_results_"+"qgwf"+"_stress_p_"+str(stress_len)+".csv"
cell_dict = {"cell_id" : cell_id_list,"reach_id":reach_ID_list}
result_df_qgwf = pd.DataFrame(cell_dict)

for time in times:
    result_qgwf = [0]*num_cells
    data = sfr.output.budget().get_data(text="GWF", totim= time)[0]["q"]
    cells_qgwf = sfr.output.budget().get_data(text="GWF", totim= time)[0]["node2"]
    
    i=0
    for j in range(0,len(result_qgwf)): 
        if cell_dict["reach_id"][j] == "No_reach":
            result_qgwf[j] = 0
            
        else:
            result_qgwf[cells_qgwf[i]-1] = data[i]
            i = i+1
        
    result_df_qgwf["q leak time "+str(time)+" (m3/s)"] = result_qgwf
result_df_qgwf.to_csv(os.path.join(result_path,csv_name))


csv_name = prjct_name+"_results_"+"stage"+"_stress_p_"+str(stress_len)+".csv"
cell_dict = {"cell_id" : cell_id_list,"reach_id":reach_ID_list}
result_df_stage = pd.DataFrame(cell_dict)

for time in times:
    i = 0
    result_stage = [0]*num_cells

    data = list(sfr.output.stage().get_data(totim = time)[0][0])
    
    for j in range(0,len(result_stage)): 
        
        if cell_dict["reach_id"][j] == "No_reach":
            result_stage[j] = 0
        else:
            result_stage[cells_qgwf[i]-1] = data[i]
            i = i+1        

    result_df_stage["h stage time "+str(time)+" (m)"] = result_stage

result_df_stage.to_csv(os.path.join(result_path,csv_name))

csv_name = prjct_name+"_results_"+"hgw"+"_stress_p_"+str(stress_len)+".csv"
cell_dict = {"cell_id" : cell_id_list,"reach_id":reach_ID_list}
result_df = pd.DataFrame(cell_dict)
i = 0
for time in times:
    result = [0]*num_cells

    result = sim.get_model().output.head().get_data(totim = time)[0][0]
    result_df["gw head time "+str(time)+" (m)"] = result

    i = i+1

result_df.to_csv(os.path.join(result_path,csv_name))









# head = hdobj.get_data()

# flows = flopy.utils.CellBudgetFile("D:/Master_Erasmus/Thesis/workSpace/mf.cbc").get_data(idx = 1)
# sfr_flows = flopy.utils.CellBudgetFile("D:/Master_Erasmus/Thesis/workSpace/mf_SFR.cbc")
# # sfr_flows.list_records()
# sfr_inpt = sfr_flows.get_data(idx = 6, text = "EXT-INFLOW")
# print(sim.get_model().output.head().get_data()[0][0])

# x = np.zeros(nvert)
# y = np.zeros(nvert)
# triangles = []

# for i in range(0,nvert):
#     x[i] = tri["vertices"][i][-2]
#     y[i] = tri["vertices"][i][-1]

# for i in range(0,ncpl):
#     triangles.append([tri["cells2d"][i][-3],tri["cells2d"][i][-2],tri["cells2d"][i][-1]])

# '''
# Simulaion plotting
# '''
# triang = mtri.Triangulation(x, y, triangles)

# fig, ax = plt.subplots()
# ax.triplot(triang, 'ko-')
# tpc = ax.tripcolor(triang,facecolors = head[0][0])
# fig.colorbar(tpc)
# plt.show()

# gwf.package_names
# gwf.get_package_list()
# u=gwf.output.methods()

#Function to get the available methods of the sfr output
# sfr_methods = sfr.output.methods()

# sfr_package_convergence = sfr.output.package_convergence()
# sfr_obs = sfr.output.obs()
# sfr_stage = sfr.output.stage()

# '''
# VTK creation
# '''
# mesh = pv.read(r"D:/Master_Erasmus/Thesis/example_basin/mesh.vtu")
# sfr_to_vtk(mesh,sfr,reach_cell_dict,name)
# simulation_result = mesh.copy()
# simulation_result.cell_data["gw_head(m)"]=sim.get_model().output.head().get_data()[0][0]
# simulation_result.save(os.path.join("D:/Master_Erasmus/Thesis/Results",'sim_gw_head.vtu'))
# fig = plt.figure(figsize=(15, 15))
# ax = plt.subplot(1, 1, 1, aspect='equal')
# img=tri.plot(ax=ax, a=head[0, 0, :], cmap='Spectral')
# fig.colorbar(img, fraction=0.02)

# fname = os.path.join(sim_ws, name + '.hds')
# hdobj = flopy.utils.HeadFile(fname, precision='double')
# head = hdobj.get_data()