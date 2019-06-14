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
import sigc4 as sigc
import tapePlacement4 as tp
import matprops_abaqus as mpa

# -----------------------------------------------------------------------------
# Define model
session.viewports['Viewport: 1'].setValues(displayedObject=None)

mdb.models.changeKey(fromName='Model-1', toName='Layer 1')
laminateModel = mdb.models['Layer 1']
laminateAssembly = laminateModel.rootAssembly

t = 2.0  # layer thickness 

# -----------------------------------------------------------------------------
# Define function used to generate parts. An object list is passed to the 
# function which then defines the geometry, the material orientations,
# the mesh, and adds the part to the assembly.
def definePart(obj, objType, objAngle, layer):
    definePart.counter += 1

    x0, y0 = obj.exterior.xy
    abqCoords = zip(x0, y0)

    profileSketch = mdb.models['Layer 1'].ConstrainedSketch(
            name='Object {}-{}'.format(layer, definePart.counter),
            sheetSize=200.0)
    for ind in range(len(abqCoords)-1):
        profileSketch.Line(
                point1=(abqCoords[ind][0], abqCoords[ind][1]),
                point2=(abqCoords[ind+1][0], abqCoords[ind+1][1]))
    
    # extrude profile to create part
    objPart = laminateModel.Part(
            name='{} {}-{}'.format(objType, layer, definePart.counter), 
            dimensionality=THREE_D, type=DEFORMABLE_BODY)
    objPart.BaseSolidExtrude(sketch=profileSketch, depth=t)
    
    # create tuple with x,y coordinates of a point within the part
    # objPoint = obj.representative_point().coords[:][0]
    objCells = objPart.cells  # all cells within part (only 1 cell)
    # define the part 'region' (necessary to define mat. orient.)
    objRegion = regionToolset.Region(cells=objCells[:])

    # TODO might need to use a set rather than a region here
    if objType == 'Resin':
        objPart.SectionAssignment(
            region=objRegion, sectionName='Resin Section',
            offset=0.0, offsetType=MIDDLE_SURFACE, offsetField='',
            thicknessAssignment=FROM_SECTION)
    else:
        objPart.SectionAssignment(
            region=objRegion, sectionName='{} Section {}'.format(objType, objAngle),
            offset=0.0, offsetType=MIDDLE_SURFACE, offsetField='',
            thicknessAssignment=FROM_SECTION)

    # mesh
    objPart.seedPart(size=1.0, deviationFactor=0.1, minSizeFactor=0.1)
    objPart.generateMesh()

    # create the part instance
    laminateAssembly = laminateModel.rootAssembly
    instanceName = '{} Instance {}-{}'.format(objType, layer, definePart.counter)
    objInstance = laminateAssembly.Instance(
            name=instanceName, part=objPart, dependent=ON)
    laminateAssembly.translate(
            instanceList=(instanceName, ), vector=(0.0, 0.0, layer*t))

    return instanceName

# ----------------------------------------------------------------------------
# Create geometries using shapely
definePart.counter = 0
laminateAngles = (0,90)
partGrid = sigc.createGrids(tapeAngles=laminateAngles, tapeWidths=(20,20))
tapePaths = tp.laminateCreation(partGrid, tapeAngles=laminateAngles, tapeWidths=(20,20), tapeSpacing=2)
partTypes = ['Tape', 'Resin', 'Undulation']

# -----------------------------------------------------------------------------
# Define materials and sections
for matAngle in laminateAngles:
    laminateModel.Material(name='Tape {}'.format(matAngle))
    rotTapeProps = mpa.rotateMatProps(140.4, 11.6, 0.289, 0.298, 6.47, 4.38, matAngle)
    laminateModel.materials['Tape {}'.format(matAngle)].Elastic(type=ENGINEERING_CONSTANTS, 
        table=(rotTapeProps, ))
    tapeSection = laminateModel.HomogeneousSolidSection(
            name='Tape Section {}'.format(matAngle), material='Tape {}'.format(matAngle))

    laminateModel.Material(name='Undulation {}'.format(matAngle))
    rotUndulationProps = mpa.rotateMatProps(112.32, 9.28, 0.289, 0.298, 5.176, 3.504, matAngle)
    laminateModel.materials['Undulation {}'.format(matAngle)].Elastic(type=ENGINEERING_CONSTANTS, 
        table=((112.32, 9.28, 9.28, 0.32, 0.32, 0.32, 5.176, 5.176, 3.504), ))
    undulationSection = laminateModel.HomogeneousSolidSection(
            name='Undulation Section {}'.format(matAngle), material='Undulation {}'.format(matAngle))

laminateModel.Material(name='Resin')
laminateModel.materials['Resin'].Elastic(type=ENGINEERING_CONSTANTS, 
        table=((7.47, 7.47, 7.47, 0.32, 0.32, 0.32, 3.5, 3.5, 3.5), ))
resinSection = laminateModel.HomogeneousSolidSection(
        name='Resin Section', material='Resin')

# -----------------------------------------------------------------------------

# create resin regions
for layerNumber in range(len(partGrid)):
    mergedObjByType = []
    objByLayer = [obj for (objID,obj) 
            in partGrid[layerNumber].iteritems()
            if obj.objectType == 'Resin']
    for objAngle in laminateAngles:
        objByAngle = [obj5 for obj5
                        in objByLayer
                        if obj5.angle == objAngle]
        # merge by angle
        mergedObjByType.append([objAngle, cascaded_union(objByAngle)])
        for obj1Angle, obj1 in mergedObjByType:           
            if obj1.geom_type == 'Polygon':
                definePart(obj1, 'Resin', obj1Angle, layerNumber)
            elif obj1.geom_type == 'MultiPolygon':
                for resinObj in obj1:
                    definePart(resinObj, 'Resin', obj1Angle, layerNumber)

# create tape regions
for (pathNumber, tapePath) in enumerate(tapePaths):
    if tapePath:
        pathAngle = tapePath[0][1].angle
        pathInstances = []
        for layerNumber in range(len(partGrid)): 
            for objType in partTypes:
                objByLayer = [obj for (objID,obj) 
                            in tapePath
                            if obj.layer == layerNumber+1]
                objByType = [obj2 for obj2
                            in objByLayer
                            if obj2.objectType == objType]
                mergedObjByType = cascaded_union(objByType)
                if mergedObjByType.geom_type == 'Polygon':
                    objInst = definePart(
                        mergedObjByType, objType, pathAngle, layerNumber)
                    pathInstances.append(objInst)
                elif mergedObjByType.geom_type == 'MultiPolygon':
                    for obj4 in mergedObjByType:
                        objInst = definePart(
                            obj4, objType, pathAngle, layerNumber)
                        pathInstances.append(objInst)
    else: continue
    if len(pathInstances) > 1:
        instanceList = [laminateAssembly.instances[name] for name in pathInstances]
        laminateAssembly.InstanceFromBooleanMerge(name='Pass {}'.format(pathNumber),
            instances=instanceList, keepIntersections=ON, originalInstances=DELETE, 
            mergeNodes=ALL, nodeMergingTolerance=1e-06, domain=BOTH)
        laminateModel.parts['Pass {}'.format(pathNumber)].seedPart(size=1.0, deviationFactor=0.1, minSizeFactor=0.1)
        laminateModel.parts['Pass {}'.format(pathNumber)].generateMesh()

laminateModel.ContactProperty('Cohesive Surface')
laminateModel.interactionProperties['Cohesive Surface'].CohesiveBehavior()
