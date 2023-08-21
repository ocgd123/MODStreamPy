# -*- coding: utf-8 -*-
"""
Created on Tue May  9 11:14:45 2023

@author: oscar
"""

def calculate_Network(river_Network, connection_dict):
    """
    Calculates river network flow and stage using provided connection information.
    
    This function calculates river network flow and stage iteratively based on the provided connection dictionary.
    
    Parameters:
        river_Network: river_Network object containing network data.
        connection_dict (dict): Dictionary containing river network connectivity information.
    
    Returns:
        river_Network: Updated river_Network object after leakage and stage calculations.
    """
    
    # Reset the river network for calculation
    river_Network.restart_Network()
    
    # Create a copy of the connection dictionary
    non_calculated_reaches = connection_dict.copy()
    # Initialize counters and lists
    orderOne_Streams_IDs = []
    total_reaches  = len(non_calculated_reaches.keys())
    solved_Reach = 0
    num_paths = 0
    solved_paths = 0
    # If no previously solved paths exist
    if len(river_Network.solved_path) == 0:
        while solved_Reach < total_reaches:
            
            solved_reaches = []
            # Iterate through reaches that haven't been calculated yet
            for dictionary_reach in non_calculated_reaches.keys():
                
                reach = dictionary_reach-1
                current_reach = river_Network.reaches[reach]
                # If the reach has only one outflow and no inflow, calculate and add flows
                if len(non_calculated_reaches[reach+1]["out"]) == 1 and len(non_calculated_reaches[reach+1 ]["in"]) == 0:
                    q_out, h_out = calculate_Reach(current_reach)
                    river_Network.reaches[reach].cv_max = calculate_Max_Conductance(current_reach,h_out)
                    river_Network.add_River_flows(reach,q_out,h_out)
                    solved_Reach = solved_Reach + 1
                    solved_reaches.append(reach)
                    continue
                # If the reach's inflow reaches are less than the number of tributaries
                if len(current_reach.in_reaches) != len(current_reach.qin): 
                    print(current_reach.reach_id)
                    continue
                
                else:
                    
                    q_out, h_out = calculate_Initial_Conditions(current_reach)
                    river_Network.reaches[reach].cv_max = calculate_Max_Conductance(current_reach,h_out)
                    river_Network.add_River_flows(reach,q_out,h_out)
                    solved_Reach = solved_Reach + 1
                    solved_reaches.append(reach)
            # If any reaches were solved in this iteration, update solved_path list
            if len(solved_reaches) > 0 :            
                for solved_reach in solved_reaches:
                    river_Network.solved_path.append(solved_reach)
                    non_calculated_reaches.pop(solved_reach+1)
    else:
        # Calculate flows and stage for previously solved paths
        for reach in river_Network.solved_path:
            current_reach = river_Network.reaches[reach]
            q_out, h_out = calculate_Reach(current_reach)
            river_Network.add_River_flows(reach,q_out,h_out)
            river_Network.calculate_RSME()
                
    return river_Network
    

def calculate_Max_Conductance(current_reach,hout):
    cv_max = current_reach.length*current_reach.width*current_reach.conductivity_aquifer/(hout-current_reach.rbot)
    return cv_max

def calculate_Initial_Conditions(reach_Section):
    """
    Calculate initial conditions for flow and stage in a reach section without leak.
    
    This function calculates the initial conditions for flow (q_down) and stage (h_down) in a reach section
    based on the provided reach section object.
    
    Parameters:
        reach_Section (reach_Section): The reach section object containing reach information.
    
    Returns:
        tuple: A tuple containing the calculated flow (q_down) and stage (h_down).
    """
        
    qin = sum(reach_Section.qin)
    if reach_Section.cross_section_geometry == "rectangle":
        
        h_down = ((qin*reach_Section.manning)/(reach_Section.width*reach_Section.slope**0.5))**(3/5)
        if h_down == 0 and reach_Section.hgw == 0:
            qgwfo = 0
            
        else:
            h_stageo = h_down + reach_Section.rbd_elev
            q_down = qin
        
        if h_down == 0 and reach_Section.hgw == 0:
            qgwf = 0
            
        else:
            h_down = ((q_down*reach_Section.manning)/(reach_Section.width*reach_Section.slope**0.5))**(3/5)
            h_stage = h_down + reach_Section.rbd_elev
        
        if reach_Section.is_in_boundary == False:
            q_down = qin
            
        if reach_Section.is_in_boundary == True:
            q_down = qin 
        
    return q_down,h_down

def calculate_Reach(reach_Section):
    """
    Calculate flow and stage in a reach section with leak.
    
    This function calculates the flow (q_down) and stage (h_down) in a reach section based on the provided reach section object.
    
    Parameters:
        reach_Section (reach_Section): The reach section object containing reach information.
    
    Returns:
        tuple: A tuple containing the calculated flow (q_down) and stage (h_down).
    """
    qin = sum(reach_Section.qin)
    if reach_Section.cross_section_geometry == "rectangle":
        
        h_down = ((qin*reach_Section.manning)/(reach_Section.width*reach_Section.slope**0.5))**(3/5)
        if h_down == 0 and reach_Section.hgw == 0:
            qgwfo = 0
            
        else:
            
            h_stageo = h_down + reach_Section.rbd_elev
            qgwfo = calculate_qL(reach_Section,h_stageo)
            
            if qgwfo > qin:
                qgwfo = qin
                
            q_down = qin - qgwfo
        
        if h_down == 0 and reach_Section.hgw == 0:
            qgwf = 0
            
        else:
            h_down = ((q_down*reach_Section.manning)/(reach_Section.width*reach_Section.slope**0.5))**(3/5)
            h_stage = h_down + reach_Section.rbd_elev
            qgwf = calculate_qL(reach_Section,h_stage)
            
            if qgwf > q_down:
                qgwf = q_down
        
        if reach_Section.is_in_boundary == False:
            reach_Section.qgwf = qgwf
            q_down = qin - reach_Section.qgwf
            
        if reach_Section.is_in_boundary == True:
            reach_Section.qgwf = 0
            q_down = qin - reach_Section.qgwf
        
    return q_down,h_down

def calculate_manning(reach_Section,q_n):
    y = ((q_n*reach_Section.manning)/(reach_Section.width*reach_Section.slope**0.5))**(3/5)
    return y

def calculateSaturation(h_riv):
    """
    Calculate saturation based on the river stage when the river is drying.
    
    This function calculates the saturation value based on the provided river stage (h_riv).
    
    Parameters:
        h_riv (float): The river stage.
    
    Returns:
        float: The calculated saturation value.
    """
    smooth = 0
    s = 1**-5
    x = h_riv
    diff = x-s
    
    if diff > 0:
        smooth = 1
        wdh = 0
        
    else: 
        aa = -1 / (s**2)
        ad = -2 / (s**2)
        b = 2 / s
        y = aa * x**2 + b * x
        dwdh = (ad * x + b)
        
        if x <= 0:
            y = 0
            dwdh = 0
            
        if diff > 1**-14:
            y = 1
            dwdh = 0
            
        smooth = y
        
    return smooth

def calculate_qL(reach_Section,h_stage):
    """
    Calculate vetical  flow exchange from between the river and the soil(qL) 
    based on river stage and groundwater depth.
    
    This function calculates the lateral groundwater flow (qL) in a reach section
    based on the provided reach section object and the specified river stage (h_stage).
    
    Parameters:
        reach_Section (reach_Section): The reach section object containing reach information.
        h_stage (float): The river stage.
    
    Returns:
        float: The calculated lateral groundwater flow (qL).
    """
    depth = h_stage-reach_Section.rbot
    sat = calculateSaturation(depth)
    # sat = 1

    if reach_Section.hgw <= reach_Section.rbot:
        qL = sat*reach_Section.conductance*(h_stage-reach_Section.rbot)

    if reach_Section.hgw > reach_Section.rbot:

        qL = sat*reach_Section.conductance*(h_stage-reach_Section.hgw) 

    return qL


def calculate_Area(reach_Sections,y):
    
    if reach_Sections.cross_section_geometry == "rectangle":
        area = reach_Sections.width * y
    