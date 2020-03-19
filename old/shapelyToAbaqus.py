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
import sigc5 as sigc
import tapePlacement5 as tp
import matprops_abaqus as mpa

# IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
# Define model name
session.viewports['Viewport: 1'].setValues(displayedObject=None)

mdb.models.changeKey(fromName='Model-1', toName='Layer 1')
laminateModel = mdb.models['Layer 1']
laminateAssembly = laminateModel.rootAssembly

# IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
# Define model geometry
# Create geometries using shapely
t = 2.0  # layer thickness 
laminateAngles = (0,90)  # define angles of tapes in laminate
partTypes = ['Tape', 'Resin', 'Undulation']
# create grid using shapely
partGrid = sigc.createGrids(tapeAngles=laminateAngles, tapeWidths=(20,20))
# identify parts in grid
tapePaths = tp.laminateCreation(
    partGrid, tapeAngles=laminateAngles, tapeWidths=(20,20), tapeSpacing=2)

# IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
# Set global parameters

# Define materials and sections
for matAngle in laminateAngles:
    laminateModel.Material(name='Tape {}'.format(matAngle))
    rotTapeProps = mpa.rotateMatProps(
        140.4, 11.6, 0.289, 0.298, 6.47, 4.38, matAngle)
    laminateModel.materials['Tape {}'.format(matAngle)].Elastic(
        type=ENGINEERING_CONSTANTS, table=(rotTapeProps, ))
    laminateModel.materials['Tape {}'.format(matAngle)].Density(
        table=((1.59e-06, ), ))
    tapeSection = laminateModel.HomogeneousSolidSection(
            name='Tape Section {}'.format(matAngle),
            material='Tape {}'.format(matAngle))
    laminateModel.Material(name='Undulation {}'.format(matAngle))
    rotUndulationProps = mpa.rotateMatProps(
        112.32, 9.28, 0.289, 0.298, 5.176, 3.504, matAngle)
    laminateModel.materials['Undulation {}'.format(matAngle)].Elastic(
        type=ENGINEERING_CONSTANTS, 
        table=((112.32, 9.28, 9.28, 0.32, 0.32, 0.32, 5.176, 5.176, 3.504), ))
    laminateModel.materials['Undulation {}'.format(matAngle)].Density(
        table=((1.59e-06, ), ))
    undulationSection = laminateModel.HomogeneousSolidSection(
            name='Undulation Section {}'.format(matAngle),
            material='Undulation {}'.format(matAngle))

laminateModel.Material(name='Resin')
laminateModel.materials['Resin'].Elastic(type=ENGINEERING_CONSTANTS, 
        table=((7.47, 7.47, 7.47, 0.32, 0.32, 0.32, 3.5, 3.5, 3.5), ))
resinSection = laminateModel.HomogeneousSolidSection(
        name='Resin Section', material='Resin')
laminateModel.materials['Resin'].Density(table=((1.1e-06, ), ))

# implicit
# create step
laminateModel.StaticStep(name='Loading Step', previous='Initial')
laminateModel.ContactProperty('Cohesive')
laminateModel.interactionProperties['Cohesive'].CohesiveBehavior(
    defaultPenalties=OFF, table=((2900.0, 2900.0, 2900.0), ))
laminateModel.interactionProperties['Cohesive'].Damage(initTable=((
    0.06, 0.06, 0.06), ), useEvolution=ON, evolutionType=ENERGY, 
    evolTable=((0.000352, ), ))


# # explicit
# # create step
# laminateModel.ExplicitDynamicsStep(name='Loading Step', previous='Initial')


# # IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
# Generate parts

# This function is used to create parts using the geometric info from Shapely
def definePart(obj, objType, objAngle, layer):
    definePart.counter += 1

    x0, y0 = obj.exterior.xy
    abqCoords = zip(x0, y0)

    profileSketch = laminateModel.ConstrainedSketch(
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

    objPart.MaterialOrientation(
        region=objRegion, orientationType=SYSTEM, axis=AXIS_2,
        localCsys=None, fieldName='', 
        additionalRotationType=ROTATION_ANGLE, additionalRotationField='',
        angle=0.0, stackDirection=STACK_3)

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
definePart.counter = 0  # initialize counter used to number parts

# ----------------------------------------------------------------------------
# In the next section the Shapely geometry objects are iterated over and
# definePart() is called to create the parts

# create resin regions
for layerNumber in range(len(partGrid)):
    resinObjByLayer = [resinObj for (m,resinObj) 
            in partGrid[layerNumber].iteritems()
            if resinObj.objectType == 'Resin']
    resinMergedByAngle = []
    for ang in laminateAngles:
        resinObjByAngle = [resinObj for resinObj 
            in resinObjByLayer
            if resinObj.angle == ang]
        resinMergedByAngle.append(cascaded_union(resinObjByAngle))
    for resinObjAngleSet in resinMergedByAngle:
        if resinObjAngleSet.geom_type == 'Polygon':
            # angle not important: resin is isotropic
            definePart(resinObjAngleSet, 'Resin', 0, layerNumber) 
        elif resinObjAngleSet.geom_type == 'MultiPolygon':
            for resinPolyObj in resinObjAngleSet:
                definePart(resinPolyObj, 'Resin', 0, layerNumber)

# create tape/undulation regions
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
                mergedObjByAngle = cascaded_union(objByType)
                if mergedObjByAngle.geom_type == 'Polygon':
                    objInst = definePart(
                        mergedObjByAngle, objType, pathAngle, layerNumber)
                    pathInstances.append(objInst)
                elif mergedObjByAngle.geom_type == 'MultiPolygon':
                    for obj4 in mergedObjByAngle:
                        objInst = definePart(
                            obj4, objType, pathAngle, layerNumber)
                        pathInstances.append(objInst)
    else: continue

    # this section merges all the parts in a tapepath into one part
    if len(pathInstances) > 1:
        instanceList = [laminateAssembly.instances[name] for name 
                in pathInstances]
        laminateAssembly.InstanceFromBooleanMerge(
            name='Pass {}'.format(pathNumber), instances=instanceList,
            keepIntersections=ON, originalInstances=DELETE, mergeNodes=ALL,
            nodeMergingTolerance=1e-06, domain=BOTH)
        laminateModel.parts['Pass {}'.format(pathNumber)].seedPart(
            size=1.0, deviationFactor=0.1, minSizeFactor=0.1)
        laminateModel.parts['Pass {}'.format(pathNumber)].generateMesh()
        passRegion = regionToolset.Region(
            cells=laminateModel.parts['Pass {}'.format(pathNumber)].cells[:])
        laminateModel.parts['Pass {}'.format(pathNumber)].MaterialOrientation(
            region=passRegion, orientationType=SYSTEM, axis=AXIS_2,
            localCsys=None, fieldName='', additionalRotationType=ROTATION_ANGLE,
            additionalRotationField='', angle=0.0, stackDirection=STACK_3)

# IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
# Create cohesive zone interactions between faces of the assembly

# implicit
laminateModel.contactDetection(defaultType=CONTACT,
        interactionProperty='Cohesive')

# # explicit
# laminateModel.ContactProperty('Tangential')
# laminateModel.interactionProperties['Tangential'].TangentialBehavior(
#     formulation=PENALTY, directionality=ISOTROPIC, slipRateDependency=OFF, 
#     pressureDependency=OFF, temperatureDependency=OFF, dependencies=0, 
#     table=((0.15, ), ), shearStressLimit=None, maximumElasticSlip=FRACTION, 
#     fraction=0.005, elasticSlipStiffness=None)
# laminateModel.ContactProperty('Cohesive')
# laminateModel.interactionProperties['Cohesive'].CohesiveBehavior()
# # laminateModel.interactionProperties['Cohesive'].Damage(initTable=((
# #     1.0, 2.0, 3.0), ), evolTable=())

# determine contacts
laminateModel.contactDetection(defaultType=CONTACT,
    interactionProperty='Cohesive', nameEachSurfaceFound=OFF,
    createUnionOfMasterSurfaces=ON, createUnionOfSlaveSurfaces=ON)
# delete interactions
for interact in laminateModel.interactions.values():
    del laminateModel.interactions['{}'.format(interact.name)]
allMaster = laminateAssembly.surfaces['CP-All-m']
allSlave = laminateAssembly.surfaces['CP-All-s']
laminateModel.ContactExp(name='GC', createStepName='Initial')
laminateModel.interactions['GC'].includedPairs.setValuesInStep(
    stepName='Initial', useAllstar=ON)
laminateModel.interactions['GC'].contactPropertyAssignments.appendInStep(
    stepName='Initial', assignments=((GLOBAL, SELF, 'Tangential'),
    (allMaster, allSlave, 'Cohesive')))

# IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
# Apply boundary conditions

# identify faces at top and bottom for coupling
tFaces = [inst.faces.getByBoundingBox(-60.0,74.0,-1.0,60.0,76.0,6.0) for inst
                in laminateAssembly.instances.values() 
                if inst.faces.getByBoundingBox(-60.0,74.0,-1.0,60.0,76.0,6.0)]
tFacesSurface = laminateAssembly.Surface(
        side1Faces=tFaces, name='Top Faces')
bFaces = [inst.faces.getByBoundingBox(-60.0,-76.0,-1.0,60.0,-74.0,6.0) for inst
                in laminateAssembly.instances.values() 
                if inst.faces.getByBoundingBox(-60.0,-76.0,-1.0,60.0,-74.0,6.0)]
bFacesSurface = laminateAssembly.Surface(
        side1Faces=bFaces, name='Bottom Faces')
symFaces = [inst.faces.getByBoundingBox(-60.0,-76.0,-1.0,60.0,76.0,0.001) for inst
                in laminateAssembly.instances.values() 
                if inst.faces.getByBoundingBox(-60.0,-76.0,-1.0,60.0,76.0,0.001)]
symFacesSet = laminateAssembly.Set(
        faces=symFaces, name='Symmetry Faces')

# create reference points for coupling                
rfPoint1Id = laminateAssembly.ReferencePoint(
        point=(0.0,80.0,2.0)).id
rfPoint1 = laminateAssembly.referencePoints[rfPoint1Id]
rfPoint1Region = laminateAssembly.Set(
        referencePoints=(rfPoint1,), name='Coupling Reference Point 1')
rfPoint2Id = laminateAssembly.ReferencePoint(
        point=(0.0,-80.0,2.0)).id
rfPoint2 = laminateAssembly.referencePoints[rfPoint2Id]
rfPoint2Region = laminateAssembly.Set(
        referencePoints=(rfPoint2,), name='Coupling Reference Point 2')

laminateModel.Coupling(
            name='Top Coupling', controlPoint=rfPoint1Region,
            surface=tFacesSurface, influenceRadius=WHOLE_SURFACE,
            couplingType=KINEMATIC, localCsys=None, u1=ON, u2=ON, u3=ON,
            ur1=ON, ur2=ON, ur3=ON)

laminateModel.Coupling(
            name='Bottom Coupling', controlPoint=rfPoint2Region,
            surface=bFacesSurface, influenceRadius=WHOLE_SURFACE,
            couplingType=KINEMATIC, localCsys=None, u1=ON, u2=ON, u3=ON,
            ur1=ON, ur2=ON, ur3=ON)

Apply boundary conditions to reference points
# implicit
laminateModel.DisplacementBC(
            name='Top Surface Constraint', createStepName='Loading Step',
            region=rfPoint1Region, u1=0.0, u2=1.0, u3=0.0, ur1=0.0, ur2=0.0,
            ur3=0.0, amplitude=UNSET, fixed=OFF, distributionType=UNIFORM,
            fieldName='', localCsys=None)
laminateModel.DisplacementBC(
            name='Bottom Surface Constraint', createStepName='Loading Step',
            region=rfPoint2Region, u1=0.0, u2=0.0, u3=0.0, ur1=0.0, ur2=0.0,
            ur3=0.0, amplitude=UNSET, fixed=OFF, distributionType=UNIFORM,
            fieldName='', localCsys=None)
# laminateModel.ZsymmBC(name='Symmetry', createStepName='Loading Step', 
#     region=symFacesSet, localCsys=None)


# # explicit
# laminateModel.SmoothStepAmplitude(name='Smoothing Amplitude',
#     timeSpan=STEP, data=((0.0, 0.0), (1e-05, 1.0)))
# laminateModel.DisplacementBC(
#     name='Bottom Surface BC', createStepName='Loading Step',
#     region=rfPoint2Region, u1=0.0, u2=0.0, u3=0.0, ur1=0.0, ur2=0.0,
#     ur3=0.0, amplitude=UNSET, fixed=OFF, distributionType=UNIFORM,
#     fieldName='', localCsys=None)
# laminateModel.VelocityBC(name='Top Surface BC', createStepName='Loading Step', 
#     region=rfPoint1Region, v1=UNSET, v2=0.01, v3=UNSET, vr1=UNSET, vr2=UNSET, 
#     vr3=UNSET, amplitude='Smoothing Amplitude', localCsys=None, 
#     distributionType=UNIFORM, fieldName='')
# laminateModel.ZsymmBC(name='Symmetry', createStepName='Loading Step', 
#     region=symFacesSet, localCsys=None)





