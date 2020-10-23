"""
This script file creates microscale models of fiber undulations.
Inputs to the model are:
- undulation amplitude h_u (mm)
- tape thickness h_f (mm)
- undulation ratio: ratio (-) 
Where the undulation ratio is the amplitude over the length of the undulation

The script creates the undulation profile in a sketch, uses this to partition
a box, splits the partitioned box into multiple parts (otherwise there can be 
no sliding between surfaces). 

Strengths and stiffnesses are given in GPa

Current model limitations:
- in the undulating region cohesive zones only defined between plies of 
  differing orientation, rather than between each ply (reasonable since
  delamination rarely occurs between plies of the same orientation)

"""

from abaqus import *
from abaqusConstants import *
import regionToolset
import mesh
import itertools
import numpy as np
import math

def rotateLaminaStiffness(C_lam, O):
    ''' This function rotates the stiffness matrix of a lamina
    by its in-plane angle (O)
    '''
    # define the transformation matrix for in-plane rotation
    Orad = math.radians(O)
    m = math.cos(Orad)
    n = math.sin(Orad)

    T_ip = np.array(
        [(m**2, n**2, 0.0, 0.0, 0.0, 2*m*n),
        (n**2, m**2, 0.0, 0.0, 0.0, -2*m*n),
        (0.0, 0.0, 1.0, 0.0, 0.0, 0.0),
        (0.0, 0.0, 0.0, m, -n, 0.0),
        (0.0, 0.0, 0.0, n, m, 0.0),
        (-m*n, m*n, 0.0, 0.0, 0.0, (m**2)-(n**2))])
    T_ip_inv = np.linalg.inv(T_ip) # invert in-plane rotation matrix

    # create the Reuter matrix
    R = np.identity(6) 
    for g in range(3,6):  # zero based indexing!
        R[g, g] = 2.0
    R_inv = np.linalg.inv(R)

    ip_rotation = np.dot(np.dot(np.dot(np.dot(
            T_ip_inv, C_lam), R), T_ip), R_inv)
    return ip_rotation

def engineeringConstants(Smatrix):
    '''This function converts a compliance matrix into a tuple of 
    engineering constants for input into a material table in Abaqus
    '''
    Ex = 1/Smatrix[0,0]
    Ey = 1/Smatrix[1,1]
    Ez = 1/Smatrix[2,2]
    nu_xy = -Smatrix[1,0]/Smatrix[0,0]
    nu_xz = -Smatrix[2,0]/Smatrix[0,0]
    nu_yz = -Smatrix[1,2]/Smatrix[1,1]
    Gxy = 1/Smatrix[5,5]
    Gxz = 1/Smatrix[4,4]
    Gyz = 1/Smatrix[3,3]
    eng_const = (Ex, Ey, Ez, nu_xy, nu_xz, nu_yz, Gxy, Gxz, Gyz)
    return eng_const

def createUndulation(h_u, h_f, ratio, theta1, theta2):
    '''
    '''

    # define initial values
    L_u = h_u/ratio  # undulation length
    L_t = 1.0 # additional refined mesh region
    w = 0.1 # specimen width
    L_g = L_u+2*L_t  # gauge length

    # initialize viewport and name model
    session.viewports['Viewport: 1'].setValues(displayedObject=None)
    mdb.Model(
            name='A-{}R-{}H-{}'.format(','.join(str(h_u).split('.')),
            ','.join(str(ratio).split('.')),','.join(str(h_f).split('.'))),
            modelType=STANDARD_EXPLICIT)
    undulationModel = mdb.models['A-{}R-{}H-{}'.format(
            ','.join(str(h_u).split('.')), ','.join(str(ratio).split('.')),
            ','.join(str(h_f).split('.')))]
    undulationAssembly = undulationModel.rootAssembly

    # <<<<<<<<<<<<<<<<<<<<<<<<<<<< CREATE PARTS >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    
    # Create the undulating fiber part
    # Sketch rectangle
    undulationProfileSketch = undulationModel.ConstrainedSketch(
            name='Undulation Profile', sheetSize=5)
    undulationProfileSketch.rectangle(point1=(-L_t, 0.0),
            point2=(L_u+L_t, h_f+h_u))

    # Create a 3D deformable part named 'Undulation' by extruding the sketch
    # the part is simply a rectangular box
    undulationPart = undulationModel.Part(
            name='Undulation', dimensionality=THREE_D, type=DEFORMABLE_BODY)
    undulationPart.BaseSolidExtrude(sketch=undulationProfileSketch, depth=w)
    undulationCells = undulationPart.cells
    undulationEdges = undulationPart.edges

    # ------------------------------------------------------------------------
    # Sketch the undulation and partition the part along its profile
    n = 10  # number of points at which the undulation is discretized
    tangent = 1.0/6.0

    # Define tangents to curve to avoid excessively thin elements
    # linear curves y = mx + b
    tan_p1 = L_u*tangent  # point at which to calculate tangent curve
    grad_p1 = ((math.pi*h_u)/(2*L_u))*np.sin((math.pi*tan_p1)/L_u) # gradient
    ytan_p1 = h_u*(np.sin((math.pi*tan_p1)/(2*L_u)))**2  # y(x) at tan_p1
    btan_p1 = ytan_p1-grad_p1*tan_p1
    tan_p1_int = (0 - btan_p1)/grad_p1  # x coord at which y = mx+b = 0.0

    tan_p2 = L_u*(1-tangent)  # point at which to calculate tangent curve
    grad_p2 = ((math.pi*h_u)/(2*L_u))*np.sin((math.pi*tan_p2)/L_u) # gradient
    ytan_p2 = h_u*(np.sin((math.pi*tan_p2)/(2*L_u)))**2  # y(x) at tan_p2
    btan_p2 = ytan_p2-grad_p2*tan_p2
    tan_p2_int = (h_u - btan_p2)/grad_p2  # x coord at which y = mx+b = h_u
            
    profile_coords = [(0.0,0.0),(tan_p1_int,0.0),(tan_p1,ytan_p1)]  
    y_sine = lambda f: h_u*(np.sin((math.pi*f)/(2*L_u)))**2
    for x in np.linspace(tan_p1,tan_p2,n):
        x_coord = x
        y_coord = y_sine(x)
        profile_coords.append((x_coord,y_coord))
    profile_coords.extend([(tan_p2,ytan_p2),(tan_p2_int,h_u),(L_u,h_u)])

    profile_upper = [(a,b+h_u) for (a,b) in profile_coords] 
    undulationProfileSketch = undulationModel.ConstrainedSketch(
            name='Undulation Profile', sheetSize=5)
    
    for k in range(len(profile_coords)-1):
        undulationProfileSketch.Line(
                point1=profile_coords[k], point2=profile_coords[k+1])
        undulationProfileSketch.Line(
                point1=profile_upper[k], point2=profile_upper[k+1])

    # connect the upper and lower profile sketches to close the geometry
    undulationProfileSketch.Line(point1=(-L_t, 0.0), point2=(-L_t, h_f))
    undulationProfileSketch.Line(point1=(L_u+L_t, h_u), point2=(L_u+L_t, h_f+h_u))
    undulationProfileSketch.Line(point1=(-L_t, 0), point2=(0, 0))
    undulationProfileSketch.Line(point1=(L_u, h_u), point2=(L_u+L_t, h_u))
    undulationProfileSketch.Line(point1=(-L_t, h_f), point2=(0, h_f))
    undulationProfileSketch.Line(point1=(L_u, h_u+h_f), point2=(L_u+L_t, h_u+h_f))

    # define the partition sketch face
    partition_sketch_face = undulationPart.faces.findAt((L_u/2.0, h_u, w), )
    # choose an edge to define the orientation of the sketch partition
    # "choose an edge vertically and to the right"
    partition_edge = undulationEdges.findAt((L_u+L_t, h_u, w), )
    # select the cell to be partitioned
    partition_cell = undulationCells.findAt((L_u/2.0, h_u, w/2.0), )

    # partition the part by loading the previous undulation profile sketch
    partition_definition = undulationPart.MakeSketchTransform(
            sketchPlane=partition_sketch_face, sketchUpEdge=partition_edge,
            sketchPlaneSide=SIDE1, origin=(0.0, 0.0, w))
    partition_sketch = undulationModel.ConstrainedSketch(
            name='Partition Sketch', sheetSize=4.54, gridSpacing=0.11,
            transform=partition_definition)
    undulationPart.projectReferencesOntoSketch(
            sketch=partition_sketch, filter=COPLANAR_EDGES)
    partition_sketch.retrieveSketch(
            sketch=undulationModel.sketches['Undulation Profile'])
    undulationPart.PartitionCellBySketch(
            sketchPlane=partition_sketch_face, sketchUpEdge=partition_edge, 
            cells=partition_cell, sketch=partition_sketch)

    # for some reason the above changes the face definition so need to create 
    # another face for effective partitioning
    partition_face = undulationPart.faces.findAt((L_u/2.0, h_u, w), )  
    # select edges defining the partition profile
    partition_boundary = ([undulationEdges[b] for b in 
            partition_face.getEdges()])
    # select edges defining the partition path (sweep path)
    partition_path = undulationEdges.findAt((L_u+L_t, h_u+h_f, w/2.0), )
    undulationPart.PartitionCellBySweepEdge(
            sweepPath=partition_path, cells=partition_cell,
            edges=partition_boundary)

    # ------------------------------------------------------------------------
    # Copy the part and then cutsweep to create two separate parts: one for the
    # undulation and one for the two sections of warp tape. This is necessary
    # because otherwise the undulation and warp surfaces can not move relative
    # to eachother

    # copy undulation part
    warpPart = undulationModel.Part(
            name='Warp', objectToCopy=undulationPart, compressFeatureList=ON)

    undulationPart.CutSweep(path=undulationEdges.findAt(((-L_t, 2*h_u, w/4.0), )),
            profile=undulationPart.faces.findAt(coordinates=(0.0, 1.5*h_u, 
                w)), flipSweepDirection=ON)

    undulationPart.CutSweep(path=undulationEdges.findAt(((-L_t, 0.0, w/4.0), )),
        profile=undulationPart.faces.findAt(coordinates=(L_u, 0.5*h_u, 
            w)), flipSweepDirection=OFF)

    # # remove undulating region from warp part
    warpPart.CutSweep(
            path=warpPart.edges.findAt(((-L_t, h_f+h_u, w/2.0), )),
            profile=warpPart.faces.findAt(
                coordinates=(L_u/2.0, (h_f+h_u)/2.0, w)),
            flipSweepDirection=ON)

    warpCells = warpPart.cells
    warpEdges = warpPart.edges

    # <<<<<<<<<<<<<<<<<<<<<<<<<<<<< MATERIALS >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    undulationMaterial = undulationModel.Material(name='Tape')
    undulationMaterial.Density(table=((1.59e-06, ), ))
    undulationMaterial.Depvar(deleteVar=20, n=24)
    undulationMaterial.UserMaterial(
            mechanicalConstants=(146.8, 11.6, 11.6, 0.289, 0.289, 0.298, 6.47, 
            6.47, 4.38, 2.61, 1.759, 0.055, 0.285, 0.105, 53.0, 0.1, 0.1, 0.00075, 
            0.0025, 0.0035))
    undulationSection = undulationModel.HomogeneousSolidSection(
            name='Tape Section', material='Tape')

    # ------------------------------------------------------------------------
    # Create sections to assign to the undulation

    all_undulation_cells = undulationPart.Set(
            cells=undulationCells[:], name='All Undulation Cells')
    all_warp_cells = warpPart.Set(cells=warpCells[:], name='All Warp Cells')

    undulationPart.SectionAssignment(
            region=all_undulation_cells, sectionName='Tape Section',
            offset=0.0, offsetType=MIDDLE_SURFACE, offsetField='',
            thicknessAssignment=FROM_SECTION)
    warpPart.SectionAssignment(
            region=all_warp_cells, sectionName='Tape Section', offset=0.0,
            offsetType=MIDDLE_SURFACE, offsetField='',
            thicknessAssignment=FROM_SECTION)

    # ------------------------------------------------------------------------
    # Assign material orientations

    # select the undulating cell and assign to region
    # undulating_cell = undulationPart.cells.findAt(((a/2.0, h_u/2.0, w/2.0), ))
    undulating_orientation_cell = regionToolset.Region(
            cells=undulationCells[:])
    warpPart.MaterialOrientation(
            region=undulating_orientation_cell, orientationType=SYSTEM, axis=AXIS_2,
            localCsys=None, fieldName='', 
            additionalRotationType=ROTATION_ANGLE, additionalRotationField='',
            angle=0.0, stackDirection=STACK_3)

    #select the two warp regions and assign material orientation
    warp_cell_region = regionToolset.Region(cells=warpCells[:])
    warpPart.MaterialOrientation(
            region=warp_cell_region, orientationType=SYSTEM, axis=AXIS_2,
            localCsys=None, fieldName='', 
            additionalRotationType=ROTATION_ANGLE, additionalRotationField='',
            angle=90.0, stackDirection=STACK_3)


if __name__ == '__main__':
    # createUndulation(h_u, h_f, ratio, theta1, theta2):
    createUndulation(0.2,0.2,0.2, 0.0, 90.0)

    L_u = L_t = 1.0
    h_f = h_u = 0.2
    w = 0.1
    undulationModel = mdb.models['A-0,2R-0,2H-0,2']
    undulationPart = undulationModel.parts['Undulation']
    undulationEdges = undulationPart.edges
