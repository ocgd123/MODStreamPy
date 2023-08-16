# -*- coding: utf-8 -*-
"""
Created on Mon Mar 27 15:10:50 2023

"""
import numpy as np
import meshio
import flopy

def vtk_to_disv(mesh,cell_type='triangle',z_interpretation='top'):
    
    """
    
    Calls a unstructured grid (e.g quads, triangles),
    
    selects all cells of a specific cell type and converts data to 2D vertices/cell
    
    information for MODFLOW, saves in one dictionary
    
     
    
    Parameters
    
    ----------
    
    file : TYPE
    
        DESCRIPTION.
    
    cell_type : TYPE, optional
    
        DESCRIPTION. The default is 'triangle'.
    
       
    
    z_interpretation: Str or None
    
        Decides how to deal with z value 'top, bottom'
    
     
    
    Returns
    
    -------
    
    disv :Dictionary
    
     
    
    """
    
    #file io
    
    if isinstance(mesh,str):
    
        print('read mesh from file')
    
        mesh=meshio.read(mesh)   
    
        #read out the cells and nodes
    
        cells=mesh.cells_dict[cell_type]
    
        points=mesh.points
    
    else:
    
        print('assume that mesh is pyvista unstructured grid')
    
        points=mesh.points
    
        if cell_type=='triangle':
    
            cells=mesh.cells_dict[5]
    
        else:
    
            raise KeyError('Only triangle mode implemented yet')
    
    disv=dict()
    
     
    
     
    
    #reconstruct cell2d : list (of lists)    innermost list contains cell number, x, y, number of vertices, and
    
        #then the vertex numbers comprising the cell.
    
     
    
    # 1. Get the centroid of each cell
    
    disv['cells2d']=list()
    
    xcycs=list()
    
    cell_z=list()
    
    for cell in range(0,cells.shape[0]):
    
        centroid=list(np.mean(points[cells[cell],0:2],axis=0))
    
        cell_z.append(np.mean(points[cells[cell],2]))
    
        cell2d=[cell]
    
        xcycs.append(centroid) # the list of centroids
    
        cell2d.extend(centroid) #centroids first
    
        cell2d.append(len(cells[cell]))
    
        cell2d.extend(list(cells[cell]))
    
        disv['cells2d'].append(cell2d)
    
    disv['xcycs']=np.array(xcycs) # to numpy array
    
    #2 get the vertices vertices : list (of lists)    innermost list contains vertex number, x, and y
    
    disv['vertices']=list()
    
       
    
    for point in range(0,mesh.points.shape[0]):
        
        disv['vertices'].append([point,points[point,0],points[point,1]])
    
    #last the number of cells per layer
    disv['ncpl']=int(cells.shape[0])
    
    #the number of vertices per layer
    disv['nvert']=int(points.shape[0])
    
    if z_interpretation is not None:
        disv[z_interpretation]=cell_z
    
    return disv

def add_mesh(self):

    """

    Returns

    -------

    None.

 

    """

    """

    if self.reload_packages and os.path.isfile(os.path.join(self.sim_dir,self.project_name+'.disv')):

        print('Reload DISV Package from simulation directory')

        return None

    else:

        self.reload_packages=False

    """
    print('add mesh to modflow')

    #convert vtk to flopy readable variables
    disv=vtk_to_disv(self.mesh,z_interpretation=self.z_interpretation)

    #how to interprete the Z Values
    if self.z_interpretation=='bottom':
        botm=disv[self.z_interpretation]
        top=(botm+np.ones((1,len(botm)))*self.aquifer_thickness).tolist()[0]

    elif self.z_interpretation=='top':
        top=disv[self.z_interpretation]
        botm=(top-np.ones((1,len(top)))*self.aquifer_thickness).tolist()[0]

    else:
        
        print('Set default by top equals thickness', self.aquifer_thickness, 'and bottom equals 0')
        top=self.aquifer_thickness
        botm=0
        
    #do the discretizsation       
    self.dis = flopy.mf6.ModflowGwfdisv(self.gwf, length_units=self.space_unit,
                                    nlay=self.nlayers, ncpl=disv['ncpl'], nvert=disv['nvert'],
                                    top=top, botm=botm,
                                    vertices=disv['vertices'], cell2d=disv['cells2d'])
    

