# -*- coding: utf-8 -*-
"""
Created on Mon May  8 09:46:48 2023

@author: oscar
"""
import warnings
import pandas as pd
import numpy as np
import os
import xml

def write_pvd (pvd_filename='test.pvd',ts_name_dict=dict()):
    """
    

    Parameters
    ----------
    pvd_filename : TYPE, optional
        DESCRIPTION. The default is 'test.pvd'.
    ts_name_dict : TYPE, optional
        DESCRIPTION. The default is dict().

    Returns
    -------
    None.

    """
        
    
    outFile = open(pvd_filename, 'w')
    
    
    pvd = xml.dom.minidom.Document()
    pvd_root = pvd.createElementNS("VTK", "VTKFile")
    pvd_root.setAttribute("type", "Collection")
    pvd_root.setAttribute("version", "0.1")
    pvd_root.setAttribute("byte_order", "LittleEndian")
    pvd.appendChild(pvd_root)
    
    collection = pvd.createElementNS("VTK", "Collection")
    pvd_root.appendChild(collection)
    
    for sim_time in ts_name_dict.keys():
        dataSet = pvd.createElementNS("VTK", "DataSet")
        dataSet.setAttribute("timestep", str(sim_time))
        dataSet.setAttribute("group", "")
        dataSet.setAttribute("part", "0")
        dataSet.setAttribute("file", ts_name_dict[sim_time])
        collection.appendChild(dataSet)
    

    pvd.writexml(outFile, newl='\n')
    outFile.close()

class reach_Section:
    """
    Represents a river reach section with various attributes and calculated variables.
    
    Parameters:
        reach_Name (str): The name of the reach.
        slope (float): The slope of the river reach.
        reach_id (int): The identifier of the reach.
        reach_cell: The cell assigned to the reach in the grid.
        manning (float): Manning's roughness coefficient for the reach.
        thickness (float): The thickness of the river section.
        length (float): The length of the river reach.
        conductivity (float): Conductivity of the river bed.
        rbd_elev (float): Elevation of the river bed.
        surface_top (float): Elevation of the groundwater surface.
        conductivity_aquifer (float): Conductivity of the aquifer.
        width (float): Width of the river section (for rectangular cross-sections).
        diversions (list, optional): List of diversion nodes if the reach has diversions.
    
    Attributes:
        is_in_boundary (bool): Whether the reach is within a head boundary.
        is_boundary (bool): Whether the reach is a boundary reach where there is inflow.
        rbot (float): Elevation of the bottom of the river bed.
        in_reaches (list): List of incoming reach sections.
        out_reaches (list): List of outgoing reach sections.
        qin (list): List of inflow values.
        hriv (float): River water depth.
        hgw (float): Groundwater head.
        qout (float): Outflow value.
        qgwf (float): Leakage flow value.
        conductance (float): Conductance of the reach.
        diversion_reaches (list): List of diversion reach sections.
        diversion_flow (list): List of diversion flow values.
    
    Methods:
        add_diversions(self, diversion): Add diversion nodes and flows to the reach.
        get_Number_Tributaries(): Get the number of tributary reaches.
        get_Number_Inflows(): Get the number of inflows to the reach.
    
    Raises:
        Exception: If required variables are missing or mismatching in diversions.
        Warning: If the reach doesn't have an assigned grid cell.
    """
        
    def __init__(self, reach_Name,slope, reach_id,
                 reach_cell,manning, 
                 thickness, length,
                 conductivity,rbd_elev,width, diversions= []):
        #Distance of the river bed to the bottom of the cell
        
        #River characteristics
        self.is_in_boundary = False
        self.rbd_elev = rbd_elev
        self.reach_name = reach_Name
        self.reach_id = reach_id
        self.reach_cell = reach_cell 
        self.slope = slope
        self.thickness = thickness
        self.length  = length
        self.manning = manning
        self.in_reaches = []
        self.out_reaches = []
        self.conductivity = conductivity
        self.rbot = rbd_elev - thickness
        self.is_boundary = False
        

        #Calculatted variables
        self.qin = []
        self.hriv = 0
        self.hgw = 0
        self.qout = 0
        self.qgwf = 0
        self.conductance = 0

        #Previous iteration variables
        self.qin_prev = []
        self.hriv_prev = 0
        self.hgw_prev = 0
        self.qout_prev = 0
        self.qgwf_prev = 0
        
        #River connecctions and flows
        self.is_Diversion  = False
        self.diversion_reaches = []
        self.diversion_flow = []
        
        #To Do: diversion control
        #To Do: Gates and control
        
        
        if slope == None or  manning == None :
            raise Exception("There is no value for a rquiered variable in reach {empty_reach}".format(empty_reach = reach_id))
            
        #Diversion data verifications
        if self.is_Diversion == True:
            
            if len(self.out_reaches) == len(diversions):
                raise Exception("Error: The diversions must not include the main reach")

            if len(self.out_reaches)-1 < len(diversions):
                raise Exception("Error: There are more diversions than reaches")                    
                
            if len(self.out_reaches)-1 > len(diversions):
                raise Exception("Error: There are missing diversions")                
                
            if len(self.diversion_flow) != len(self.diversion_nodes):
                raise Exception("There is no value for a rquiered variable in reach {empty_reach}".format(empty_reach = reach_id))
            
        if len(diversions)>0:
            
            self.is_Diversion  = True
            self.diversion_nodes = diversions
            

            
        if reach_cell == None:
            warnings.warn("Warning: The reach {empty_reach} has no cell assigned in the grid".format(empty_reach = reach_id))
        
        self.width = width
        self.conductance = (self.width*self.conductivity*self.length)/self.thickness
            
            
        def add_diversions(self,diversion):
            """
            Add diversion nodes and flows to the reach.
            
            Parameters:
                diversion (list): List containing diversion information.
            
            Returns:
                None
            """
                        
            self.is_Diversion = True
            diversion[0]
            
        def get_Number_Tributaries():
            """
            Get the number of tributary reaches.
            
            Returns:
                int: Number of tributary reaches.
            """
            return len(self.tributaries)
        
        def get_Number_Inflows():
            """
            Get the number of inflows to the reach.
    
            Returns:
                int: Number of inflows to the reach.
            """
            return len(self.qin)
        
            
class river_Network:
    
    """
    Represents a network of river reaches with attributes for managing flows and calculations.
    
    Attributes:
        reaches (dict): Dictionary containing reach identifiers as keys and reach_Section objects as values.
        solved_path (list): List of solved path elements in the network.
        rsme_hriv (float): Root Mean Squared Error for river stage height.
        rsme_hgw (float): Root Mean Squared Error for groundwater head.
        rsme_qout (float): Root Mean Squared Error for outflow values.
        rsme_qgwf (float): Root Mean Squared Error for groundwater flow values.
        qin_total (float): Total inflow from boundaries.
        qout_total (float): Total outflow to boundaries.
        leak_total (float): Total groundwater leakage.
    
    Methods:
        create_Network(self, ...): Creates a river network based on provided data.
        add_River_flows(self, reach_id, q, h): Adds river flows to the network.
        restart_Network(self):         Restarts the network by storing previous iteration variables and setting
        the inlet of the reaches that are no boundaries to 0.
        update_gw_Heads(self, sim, modflow_time): Updates groundwater heads in the network.
        calculate_RSME(self): Calculates Root Mean Squared Errors for various variables.
        calculate_Water_Balance(self): Calculates total inflow, leakge and outflow for water balance.
        add_Boundaries(self, inflow_dict): Adds inflow boundaries to the network.
        add_Reach_in_Boundaries(self, boundary_cells): Flags reaches within boundary cells.
    
    """
    
    
    def __init__(self):
        self.reaches = {}
        self.solved_path = []
        self.rsme_hriv = 100
        self.rsme_hgw = 100
        self.rsme_qout = 100
        self.rsme_qgwf = 100
        self.qin_total = 0
        self.qout_total = 0
        self.leak_total = 0

        
    def create_Network(self, name_list,slope_list,reach_id_list,
                       reach_cell_list,
                       conect_data, length_list, 
                       rbd_elev_list,
                       width_list,riv_hydr_prop,
                       diversion_df = pd.DataFrame()):
        """
        Creates a river network based on provided data.
        
        Parameters:
            name_list (list): List of reach names.
            slope_list (list): List of reach slopes.
            reach_id_list (list): List of reach identifiers.
            ... (other parameters)
        
        Returns:
            None
        """
        
        for i in range(0,len(name_list)):
            
            manning = riv_hydr_prop[name_list[i]]["River properties"]["river_manning"]
            conductivity = riv_hydr_prop[name_list[i]]["River properties"]["river_conductivity"]
            thick = riv_hydr_prop[name_list[i]]["River properties"]["river_bed_thickness"]
            
            reach_section = reach_Section(name_list[i], slope_list[i], 
                                          reach_id_list[i], reach_cell_list[i],
                                          manning,thick,
                                          length_list[i],conductivity,
                                          rbd_elev_list[i], 
                                          width_list[i])
            
            reach_section.in_reaches = conect_data[reach_id_list[i]+1]["in"]
            reach_section.out_reaches = [abs(out_id) for out_id in conect_data[reach_id_list[i]+1]["out"]]
            
            self.reaches[reach_id_list[i]] = reach_section
        
        if len(diversion_df)>0:
            
            for i in range(0,len(diversion_df)):
                reach = diversion_df.loc[i][1]-1
                div_reach = diversion_df.loc[i][2]
                div_type = diversion_df.loc[i][3]
                q = diversion_df.loc[i][4]
                self.reaches[reach].is_Diversion = True
                
                if div_type == "diversion":
                    self.reaches[reach].diversion_reaches.append(div_reach-1)
                    self.reaches[reach].diversion_flow.append(q)



    def add_River_flows(self, reach_id, q, h):
        """
        Adds river flows to the network.
        
        Parameters:
            reach_id (int): Identifier of the reach.
            q (float): Flow value.
            h (float): River stage height.
        
        Returns:
            None
        """
        if self.reaches[reach_id].is_Diversion == False:
            
            self.reaches[reach_id].qout = q
            self.reaches[reach_id].hriv = h
            out_reaches = self.reaches[reach_id].out_reaches
            if out_reaches != 0:
                for out_reach in out_reaches:
                    self.reaches[out_reach].qin.append(q)
                        
                        
        if self.reaches[reach_id].is_Diversion == True:
            
            self.reaches[reach_id].qout = q
            self.reaches[reach_id].hriv = h
            out_reaches = self.reaches[reach_id].out_reaches
            div_flows = self.reaches[reach_id].diversion_flow
            div_reaches = self.reaches[reach_id].diversion_reaches

            if out_reaches != 0:
                for out_reach in out_reaches:
                    if out_reach not in div_reaches:
                        self.reaches[out_reach].qin.append(q-sum(div_flows))
                        break
                    
            for i in range(0, len(div_flows)):
                out_reach = div_reaches[i]
                out_flow = div_flows[i]
                self.reaches[out_reach].qin.append(out_flow)
            
    def restart_Network(self):
        """
        Restarts the network by storing previous iteration variables and setting
        the inlet of the reaches that are no boundaries to 0.
          
        Returns:
            None
        """  
            
        for reach in self.reaches.keys():
            
            self.reaches[reach].qin_prev = self.reaches[reach].qin
            self.reaches[reach].hriv_prev = self.reaches[reach].hriv
            self.reaches[reach].qout_prev = self.reaches[reach].qout
            self.reaches[reach].qgwf_prev = self.reaches[reach].qgwf
            
            if self.reaches[reach].is_boundary == False: 
                self.reaches[reach].qin = []


    def update_gw_Heads(self,sim,modflow_time):
        """
        Updates groundwater heads in the network.
        
        Parameters:
            sim: Simulation object.
            modflow_time: Time value for MODFLOW output.
        
        Returns:
            None
        """
        head = sim.get_model().output.head().get_data(totim = modflow_time)[0][0]
        for reach in self.reaches.keys():
            gw_cell = self.reaches[reach].reach_cell
            self.reaches[reach].hgw_prev = self.reaches[reach].hgw
            self.reaches[reach].hgw = head[gw_cell]    
            
    def calculate_RSME(self):
        """
        Calculates Root Mean Squared Errors for various variables.
        
        Returns:
            None
        """
        rsme_hriv = []
        rsme_hgw =[]
        rsme_qout = []
        rsme_qgwf = []
        
        for reach in self.reaches.keys():
            rsme_hriv.append((self.reaches[reach].hriv_prev - self.reaches[reach].hriv)**2)
            rsme_hgw.append((self.reaches[reach].hgw_prev - self.reaches[reach].hgw)**2)
            rsme_qout.append((self.reaches[reach].qout_prev - self.reaches[reach].qout)**2)
            rsme_qgwf.append((self.reaches[reach].qgwf_prev - self.reaches[reach].qgwf)**2)
            
        self.rsme_hriv = sum(rsme_hriv)/len(rsme_hriv)
        self.rsme_hgw = sum(rsme_hgw)/len(rsme_hgw)
        self.rsme_qout = sum(rsme_qout)/len(rsme_qout)
        self.rsme_qgwf = sum(rsme_qgwf)/len(rsme_qgwf)
        
    def calculate_Water_Balance(self):
        """
        Calculates total inflow and outflow for water balance.
        
        Returns:
            None
        """
        qgwf =[]
        qout = []
        for reach in self.reaches.keys():
            self.leak_total = self.reaches[reach].qgwf+ self.leak_total 
            if len(self.reaches[reach].out_reaches) == 0:
                self.qout_total = self.qout_total + self.reaches[reach].qout
     
    def add_Boundaries(self, inflow_dict):
        """
        Adds inflow boundaries to the network.
        
        Parameters:
            inflow_dict (dict): Dictionary of inflow values keyed by reach identifiers.
        
        Returns:
            None
        """
        for inflow_ID in inflow_dict:
            self.qin_total = self.qin_total + inflow_dict[inflow_ID]
            self.reaches[inflow_ID].qin.append(inflow_dict[inflow_ID])
            self.reaches[inflow_ID].is_boundary = True
            
    def add_Reach_in_Boundaries(self, boundary_cells):
        """
        Flags reaches within boundary cells.
        
        Parameters:
            boundary_cells (list): List of boundary cell IDs.
        
        Returns:
            None
        """
        
        for reach in self.reaches.keys():
            reach_cell = self.reaches[reach].reach_cell
            if reach_cell in boundary_cells:
                self.reaches[reach].is_in_boundary = True

class time_Manager:
    """
    Manages time-related information for a MODFLOW simulation.
    
    Parameters:
        sim: MODFLOW simulation object.
    
    Attributes:
        modflow_times (list): List of lists containing MODFLOW stress period times.
        modflow_time_steps (list): List of lists containing MODFLOW stress period time steps.
        modflow_stperiod_len (list): List of MODFLOW stress period lengths.
        modflow_times_seconds (list): List of lists containing stress period times converted to seconds.
        modflow_time_steps_seconds (list): List of lists containing stress period time steps converted to seconds.
        gw_units (str): Groundwater time units from the simulation.
        modflow_stperiods (list): List of MODFLOW stress period indices.
    
    Methods:
        convert_Time_to_Seconds(self, unit, value): Converts a time value from specified unit to seconds.
        convert_from_Seconds(self, unit, value): Converts a time value from seconds to specified unit.
    
    """
    def __init__(self, sim):
        """
        Initializes the time_Manager object.

        Parameters:
            sim: MODFLOW simulation object.

        Returns:
            None
        """
        modflow_time_info = sim.get_package("tdis").blocks["perioddata"].datasets['perioddata'].get_data()
        self.modflow_times = []
        self.modflow_time_steps = []
        self.modflow_stperiod_len = []
        
        self.modflow_times_seconds = []
        self.modflow_time_steps_seconds = []
        self.modflow_stperiod_len_seconds = []
        
        self.gw_units = sim.get_package("tdis").blocks["options"].datasets["time_units"].data
        stress_time = 0
        for modflow_stperiod in modflow_time_info:
            
            stressp_len= modflow_stperiod[0]
            stressp_time_interval= modflow_stperiod[1]
            stressp_mult= modflow_stperiod[2]
            
            if stressp_mult == 1:
                dt1 = stressp_len/stressp_time_interval

            if stressp_mult != 1:
                dt1 = stressp_len*(stressp_mult-1)/(stressp_mult**stressp_time_interval-1)
            
            dtold =dt1
            tacum = dt1
            stress_time_steps = []
            stress_times = []
            self.modflow_stperiod_len.append(stressp_len)
            stress_time_steps.append(dt1)
            stress_times.append(tacum+stress_time)
            
            while tacum != stressp_len:
                dtnew = dtold*stressp_mult
                tacum = tacum+dtnew
                stress_time_steps.append(dtnew)
                stress_times.append(tacum+stress_time)
                dtold = dtnew
                
            self.modflow_time_steps.append(stress_time_steps)
            self.modflow_times.append(stress_times)
            stress_time = stress_time+stressp_len

        self.modflow_stperiods = list(range(0,len(self.modflow_stperiod_len)))
        
        for stperiod in self.modflow_stperiod_len:
            stperiod_len_sec = self.convert_Time_to_Seconds(self.gw_units, stperiod)
            self.modflow_stperiod_len_seconds.append(stperiod_len_sec)
            
        for stperiod in self.modflow_stperiods:
            
            time_steps_sec = []
            for time_step in self.modflow_time_steps[stperiod]:
                time_step_sec = self.convert_Time_to_Seconds(self.gw_units, time_step)
                time_steps_sec.append(time_step_sec)
            self.modflow_time_steps_seconds.append(time_steps_sec)
                
            times_sec = []
            for time in self.modflow_times[stperiod]:
                time_sec = self.convert_Time_to_Seconds(self.gw_units, time)
                times_sec.append(time_sec)
            self.modflow_times_seconds.append(times_sec)
            
    def convert_Time_to_Seconds(self,unit, value):
        """
        Converts a time value from specified unit to seconds.
        
        Parameters:
            unit (str): Time unit to be converted from.
            value (float): Time value to be converted.
        
        Returns:
            float: Time value converted to seconds.
        """
        if  unit == "years":
            value_converted = value*365*24*60*60
            
        elif unit== "days":
            value_converted = value*24*60*60
            
        elif unit== "hours":
            value_converted = value*60*60
            
        elif unit == "minutes":
            value_converted = value*60
            
        elif unit == "seconds":           
            value_converted = value
            
        else:
            raise TypeError("Please introduce a valid unit")
        
        return value_converted
    
    def convert_from_Seconds(self, unit,value):
        """
        Converts a time value from seconds to specified unit.
        
        Parameters:
            unit (str): Time unit to be converted to.
            value (float): Time value to be converted.
        
        Returns:
            float: Time value converted from seconds to specified unit.
        """
        if  unit == "years":
            value_converted = value/(365*24*60*60)
            
        elif unit== "days":
            value_converted = value/(24*60*60)
            
        elif unit== "hours":
            value_converted = value/(60*60)
            
        elif unit == "minutes":
            value_converted = value/(60)
            
        elif unit == "seconds":           
            value_converted = value
            
        else:
            raise TypeError("Please introduce a valid unit")
        return value_converted
            


class results_Manager:
    """
    Manages and exports simulation results.
    
    Parameters:
        modflow_stperiods (list): List of MODFLOW stress period indices.
        modflow_times (list): List of lists containing MODFLOW stress period times.
        num_cells (int): Number of cells in the simulation.
        result_path (str): Path to save the results.
    
    Attributes:
        qin_noleak (list): List of inflow values without leakage.
        hriv_noleak (list): List of river stage values without leakage.
        hgw_noleak (list): List of groundwater head values without leakage.
        qout_noleak (list): List of outflow values without leakage.
        qgwf_noleak (list): List of groundwater flow values without leakage.
        water_balance (dict): Dictionary tracking water balance for each stress period.
        num_cells (int): Number of cells in the simulation.
        results_dict (dict): Dictionary containing simulation results.
        results_dict_prev (dict): Dictionary containing previous iteration simulation results.
        result_path (str): Path to save the results.
    
    Methods:
        add_no_Leak_Network(self, river_Network): Adds network results without leakage.
        export_Data_csv(self, prjct_name, river_Network): Exports simulation results as CSV files.
        export_Data_vtu(self, mesh, prjct_name): Exports simulation results as VTK (VTU) files.
        add_Initial_Stage(self, river_stage_initial): Adds initial river stage data.
        update_qLeak(self, Stress_Periods, qleak_list): Updates leakage flow data.
        update_head(self, Stress_Periods, gw_head_list): Updates groundwater head data.
        update_stage(self, Stress_Periods, river_stage_list): Updates river stage data.
    """
    def __init__(self,modflow_stperiods,modflow_times, num_cells,result_path): 
        """
        Initializes the results_Manager object.
        
        Parameters:
            modflow_stperiods (list): List of MODFLOW stress period indices.
            modflow_times (list): List of lists containing MODFLOW stress period times.
            num_cells (int): Number of cells in the simulation.
            result_path (str): Path to save the results.
        
        Returns:
            None
        """
        
        self.qin_noleak = []
        self.hriv_noleak = []
        self.hgw_noleak = []
        self.qout_noleak = []
        self.qgwf_noleak = []
        self.water_balance = {}
        
        self.num_cells  =num_cells
        self.results_dict = {}
        self.results_dict_prev = {}

        self.result_path = result_path
        for stress in modflow_stperiods:
            
            self.water_balance[stress] = 0
            self.results_dict[stress] = {}
            self.results_dict_prev[stress] = {}
            self.water_balance[stress] =  0
            
            for time in modflow_times:
                
                qin = dict((period,[]) for period in time)
                hriv = dict((period,[]) for period in time)
                hgw = dict((period,[]) for period in time)
                qout = dict((period,[]) for period in time)
                qgwf = dict((period,[]) for period in time)
                
                qin_prev = dict((period,[]) for period in time)
                hriv_prev = dict((period,[]) for period in time)
                hgw_prev = dict((period,[]) for period in time)
                qout_prev = dict((period,[]) for period in time)
                qgwf_prev = dict((period,[]) for period in time)
                
            self.results_dict[stress]["qin"] = qin
            self.results_dict[stress]["hriv"] = hriv
            self.results_dict[stress]["hgw"] = hgw
            self.results_dict[stress]["qout"] = qout
            self.results_dict[stress]["qgwf"] = qgwf
            
            self.results_dict_prev[stress]["qin"] = qin_prev
            self.results_dict_prev[stress]["hriv"] = hriv_prev
            self.results_dict_prev[stress]["hgw"] = hgw_prev
            self.results_dict_prev[stress]["qout"] = qout_prev
            self.results_dict_prev[stress]["qgwf"] = qgwf_prev    

            
    def add_no_Leak_Network(self, river_Network):
        """
        Adds network results without leakage.
        
        Parameters:
            river_Network: river_Network object containing network data.
        
        Returns:
            None
        """
        qin_noleak = np.zeros(self.num_cells)
        hriv_noleak = np.zeros(self.num_cells)
        qout_noleak = np.zeros(self.num_cells)
        for reach in river_Network.reaches.keys():
            qin_noleak[river_Network.reaches[reach].reach_cell] = sum(river_Network.reaches[reach].qin)
            hriv_noleak[river_Network.reaches[reach].reach_cell] = river_Network.reaches[reach].hriv
            qout_noleak[river_Network.reaches[reach].reach_cell] = river_Network.reaches[reach].qout

        self.qin_noleak = qin_noleak
        self.hriv_noleak = hriv_noleak
        self.qout_noleak = qout_noleak
    
    def export_Data_csv(self,prjct_name,river_Network):
        """
        Exports simulation results as CSV files.
        
        Parameters:
            prjct_name (str): Name of the project.
            river_Network: river_Network object containing network data.
        
        Returns:
            None
        """
        reach_ID_list = ["No_reach"] * self.num_cells
        for reach in river_Network.reaches.keys():
            reach_ID_list[river_Network.reaches[reach].reach_cell] = river_Network.reaches[reach].reach_id+1
        
        cell_id = []
            
            
        for stress_data in self.results_dict:
            data_types = self.results_dict[stress_data]
            
            for data_type in data_types:
                times = data_types[data_type]
                cell_dict = {"cell_id" : range(1, self.num_cells+1),"reach_id":reach_ID_list}
                result_df = pd.DataFrame(cell_dict)
        
                csv_name = prjct_name+"_results_"+data_type+"_stress_p"+str(stress_data)+".csv"
                for time in times:
                    
                    if data_type == "hriv":
                        data_label = "h stage time "+str(time)+" (m)"
                        data = self.results_dict[stress_data][data_type][time]
                        result_df[data_label] = data
                        
                    if data_type == "qgwf":
                        data_label = "q leak time "+str(time)+" (m3/s)"
                        data = self.results_dict[stress_data][data_type][time] 
                        result_df[data_label] = data

                    if data_type == "hgw":
                        data_label = "gw head time "+str(time)+" (m)"
                        data = self.results_dict[stress_data][data_type][time] 
                        result_df[data_label] = data
                        
                    if data_type == "qin":
                        data_label = "reach initial node qin"+str(time)+" (m3/s)"
                        data = self.results_dict[stress_data][data_type][time] 
                        result_df[data_label] = data
                    
                    if data_type == "qout":
                        data_label = "reach final node qout"+str(time)+" (m3/s)"
                        data = self.results_dict[stress_data][data_type][time] 
                        result_df[data_label] = data
                        
                result_df.to_csv(os.path.join(self.result_path,csv_name))

    def export_Data_vtu(self,mesh,prjct_name):
        """
        Exports simulation results as VTK (VTU) files.
        
        Parameters:
            mesh: Mesh object containing simulation mesh data.
            prjct_name (str): Name of the project.
        
        Returns:
            None
        """
        all_results = mesh.copy()
        for stress_data in self.results_dict:   
            data_types = self.results_dict[stress_data]
            
            for data_type in data_types:
                times = data_types[data_type]
                vtk_names = {}
                for time in times:

                    all_results = mesh.copy()
                    vtk_name = prjct_name+"_"+data_type+"_"+str(time)+".vtu"
                    vtk_names[time] = vtk_name
                    
                    if data_type == "hriv":
                        all_results.cell_data["h_stage_(m)"]=self.results_dict[stress_data][data_type][time]
                        all_results.save(os.path.join(self.result_path,vtk_name))
                        
                    if data_type == "qgwf":

                        all_results.cell_data["q_leak_(m3/s)"]=self.results_dict[stress_data][data_type][time]
                        all_results.save(os.path.join(self.result_path,vtk_name))
                        
                    if data_type == "hgw":
                        all_results.cell_data["gw_head_(m)"] = self.results_dict[stress_data][data_type][time]
                        all_results.save(os.path.join(self.result_path,vtk_name))
                        
                pvd_name = prjct_name+"_"+data_type+"_stress_p_"+str(stress_data)+".pvd"
                path = os.path.join(self.result_path,pvd_name)
                write_pvd(path,vtk_names)
            # simulation_result = mesh.copy()
            # simulation_result.cell_data["gw_head (m)"]=head
            

            
    def add_Initial_Stage(self, river_stage_initial):
        """
        Adds initial river stage data.
        
        Parameters:
            river_stage_initial: Initial river stage data.
        
        Returns:
            None
        """
        self.river_stage_initial  = {}
        
    def update_qLeak(self, Stress_Periods, qleak_list):
        """
        Updates leakage flow data.
        
        Parameters:
            Stress_Periods: Stress period index.
            qleak_list (list): Leakage flow data.
        
        Returns:
            None
        """
        self.qleak[Stress_Periods] = qleak_list
        
    def update_head(self, Stress_Periods, gw_head_list):
        """
        Updates groundwater head data.
        
        Parameters:
            Stress_Periods: Stress period index.
            gw_head_list (list): Groundwater head data.
        
        Returns:
            None
        """
        self.gw_head[Stress_Periods] = gw_head_list
        
    def update_stage(self, Stress_Periods, river_stage_list):
        """
        Updates river stage data.
        
        Parameters:
            Stress_Periods: Stress period index.
            river_stage_list (list): River stage data.
        
        Returns:
            None
        """
        self.river_stage[Stress_Periods] = river_stage_list