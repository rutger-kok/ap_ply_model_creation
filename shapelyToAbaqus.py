from sys import path
path.append(r'C:\Python27\Lib\site-packages')
path.append(r"\\arran.sms.ed.ac.uk\home\s1342398\GitHub\interlaced_model_creation")   
from abaqus import *
from abaqusConstants import *
import mesh
import regionToolset
import math
import job
import displayGroupMdbToolset as dgm
from shapely.geometry import *
from shapely import affinity
import numpy as np
from shapely.ops import cascaded_union
from itertools import count
from shapelyInterlacedGeometryCreation import *

# -----------------------------------------------------------------------------
# Define model
session.viewports['Viewport: 1'].setValues(displayedObject=None)

mdb.models.changeKey(fromName='Model-1', toName='Layer 1')
laminateModel = mdb.models['Layer 1']
laminateAssembly = laminateModel.rootAssembly

t = 2.0  # layer thickness 
# -----------------------------------------------------------------------------
# Define materials and sections
laminateModel.Material(name='Tape')
laminateModel.materials['Tape'].Elastic(type=ENGINEERING_CONSTANTS, 
        table=((140.4, 11.6, 11.6, 0.289, 0.289, 0.298, 6.47, 6.47, 4.38), ))
laminateModel.Material(name='Resin')
laminateModel.materials['Resin'].Elastic(type=ENGINEERING_CONSTANTS, 
        table=((7.47, 7.47, 7.47, 0.32, 0.32, 0.32, 3.5, 3.5, 3.5), ))


tapeSection = laminateModel.HomogeneousSolidSection(
        name='Tape Section', material='Tape')
resinSection = laminateModel.HomogeneousSolidSection(
        name='Resin Section', material='Resin')

# -----------------------------------------------------------------------------
# Define function used to generate parts. An object list is passed to the 
# function which then defines the geometry, the material orientations,
# the mesh, and adds the part to the assembly.
def definePart(objectList, objType):
    for num, obj in enumerate(objectList):
        if obj.area > 0.0001:
            # Identify object type
            # pass

            # Extract coordinates and sketch profile
            x0, y0 = obj.exterior.xy
            x0r = map(lambda x: round(x, 5), x0) 
            y0r = map(lambda x: round(x, 5), y0)
            abqCoords = zip(x0r, y0r)

            profileSketch = mdb.models['Layer 1'].ConstrainedSketch(
                    name='Object {}'.format(num), sheetSize=200.0)
            for ind in range(len(abqCoords)-1):
                profileSketch.Line(
                        point1=(abqCoords[ind][0], abqCoords[ind][1]),
                        point2=(abqCoords[ind+1][0], abqCoords[ind+1][1]))
            
            # Extrude profile to create part
            objPart = laminateModel.Part(
                    name='{} {}'.format(objType, num), dimensionality=THREE_D, 
                    type=DEFORMABLE_BODY)
            objPart.BaseSolidExtrude(sketch=profileSketch, depth=t)
            
            
            # create tuple with x,y coordinates of a point within the part
            objPoint = obj.representative_point().coords[:][0]
            objCells = objPart.cells  # all cells within part (only 1 cell)
            # define the part 'region' (necessary to define mat. orient.)
            objRegion = regionToolset.Region(cells=objCells[:])

            # TODO might need to use a set rather than a region here
            objPart.SectionAssignment(
                    region=objRegion, sectionName='{} Section'.format(objType),
                    offset=0.0, offsetType=MIDDLE_SURFACE, offsetField='',
                    thicknessAssignment=FROM_SECTION)

            # rotate material orientations to coincide with Tape angle
            objPart.MaterialOrientation(
                    region=objRegion, orientationType=SYSTEM, axis=AXIS_3,
                    localCsys=None, fieldName='', 
                    additionalRotationType=ROTATION_ANGLE,
                    additionalRotationField='', angle=obj.angleLabel,
                    stackDirection=STACK_3)

            # mesh
            # create instance


            # define INTERACTIONS! might be an issue

        else: continue

# define step

# define BCs


test = definePart(list(Tape.getinstances(1)), 'Tape')
test2 = definePart(list(Resin.getinstances(1)), 'Resin')

# from matprops_abaqus import *
# assign material properties
# rotatedMaterialProps = matprops(E_11=140.4,E_22=11.6,nu_12=0.289,
#         nu_23=0.298,G_12=6.47,G_23=4.38, angle=obj.angleLabel)

# objTStr = str(type(obj))  # object type string
# objType = objTStr[objTStr.rfind('.')+1: objTStr.rfind('\'')]