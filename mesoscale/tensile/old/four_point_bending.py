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
import chou1972v2 as chou
import sigc10 as sigc
import tapePlacement10 as tp
import matprops_abaqus as mpa

# IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
# Define model name
session.viewports['Viewport: 1'].setValues(displayedObject=None)

mdb.models.changeKey(fromName='Model-1', toName='Layer 1')
laminateModel = mdb.models['Layer 1']
laminateAssembly = laminateModel.rootAssembly

# IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
# Define model geometry

cpt = 0.175  # cured ply thickness 
laminateAngles = (0,90)  # define angles of tapes in laminate
tw = (25,25)  # tape widths
# Define undulation unit cell geometry
n = 200  # number of increments
ht = 0.35  # total height of the laminate
hf = 0.175  # height of undulating plies (ply thickness * num of plies)
hu = 0.175  # height of the undulation
uw = 2.5  # total length of the undulation
O1 = math.radians(laminateAngles[1])  # in-plane angle of non-undulating plies
O2 = math.radians(laminateAngles[0])  # in-plane angle of undulating plies
tLaminate = cpt*len(laminateAngles)

partTypes = ['Tape', 'Resin', 'Undulation']
# create grid using shapely
partGrid = sigc.createGrids(tapeAngles=laminateAngles, tapeWidths=tw,
                            undulationWidth=uw)
# identify parts in grid
tapePaths = tp.laminateCreation(
    grid=partGrid, tapeAngles=laminateAngles, tapeWidths=tw,
    tapeSpacing=1, undulationWidth=uw)

# IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
# Set global parameters

Clamina = chou.CFromConstants(140.4,11.6,0.289,0.298,6.47,4.38)
Vfrac, a, rat = chou.undulationGeometry(ht, hf, hu, uw, n)
CIsostrain, CIsostress = chou.determineStiffness(Clamina, a, O1, O2, Vfrac, n)

laminateModel.Material(name='Tape')
rotTapeProps = mpa.rotateMatProps(140.4, 11.6, 0.289, 0.298, 6.47, 4.38, 0.0)
laminateModel.materials['Tape'].Elastic(
    type=ENGINEERING_CONSTANTS, table=(rotTapeProps, ))
laminateModel.materials['Tape'].Density(table=((1.59e-06, ), ))
tapeSection = laminateModel.HomogeneousSolidSection(
        name='Tape Section', material='Tape')

laminateModel.Material(name='Undulation')
rotUndulationProps = chou.engineeringConstants(
    np.linalg.inv(chou.rotateLaminaStiffness(CIsostress, math.radians(0.0), 0)))
laminateModel.materials['Undulation'].Elastic(
    type=ENGINEERING_CONSTANTS, table=(rotUndulationProps, ))
laminateModel.materials['Undulation'].Density(table=((1.59e-06, ), ))
undulationSection = laminateModel.HomogeneousSolidSection(
        name='Undulation Section', material='Undulation')

laminateModel.Material(name='Resin')
laminateModel.materials['Resin'].Elastic(type=ENGINEERING_CONSTANTS, 
        table=((7.47, 7.47, 7.47, 0.32, 0.32, 0.32, 3.5, 3.5, 3.5), ))
resinSection = laminateModel.HomogeneousSolidSection(
        name='Resin Section', material='Resin')
laminateModel.materials['Resin'].Density(table=((1.1e-06, ), ))

# explicit
# create step
laminateModel.ExplicitDynamicsStep(name='Loading Step', previous='Initial',
        timePeriod=30.0)
laminateModel.fieldOutputRequests['F-Output-1'].setValues(
        numIntervals=100)

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
    objPart.BaseSolidExtrude(sketch=profileSketch, depth=cpt)
    
    # create tuple with x,y coordinates of a point within the part
    # objPoint = obj.representative_point().coords[:][0]
    objCells = objPart.cells  # all cells within part (only 1 cell)
    # define the part 'region' (necessary to define mat. orient.)
    objRegion = regionToolset.Region(cells=objCells[:])

    if objType == 'Resin' or objType == 'Tape':
        objPart.SectionAssignment(
            region=objRegion, sectionName='Tape Section',
            offset=0.0, offsetType=MIDDLE_SURFACE, offsetField='',
            thicknessAssignment=FROM_SECTION)
    else:
        objPart.SectionAssignment(
            region=objRegion, sectionName='Undulation Section',
            offset=0.0, offsetType=MIDDLE_SURFACE, offsetField='',
            thicknessAssignment=FROM_SECTION)

    objPart.MaterialOrientation(
        region=objRegion, orientationType=SYSTEM, axis=AXIS_3,
        localCsys=None, fieldName='', 
        additionalRotationType=ROTATION_ANGLE, additionalRotationField='',
        angle=objAngle, stackDirection=STACK_3)

    # mesh
    objPart.seedPart(size=1.0, deviationFactor=0.1, minSizeFactor=0.1)
    objPart.generateMesh()

    # mirror part
    mirrorPart = laminateModel.Part(
            name='{} {}-{} Mirror'.format(objType, layer, definePart.counter), 
            objectToCopy=laminateModel.parts['{} {}-{}'.format(
                objType, layer, definePart.counter)],
            compressFeatureList=ON, mirrorPlane=XYPLANE)

    mirrorCells = mirrorPart.cells  # all cells within part (only 1 cell)
    # define the part 'region' (necessary to define mat. orient.)
    mirrorRegion = regionToolset.Region(cells=mirrorCells[:])

    if objType == 'Resin' or objType == 'Tape':
        mirrorPart.SectionAssignment(
            region=mirrorRegion, sectionName='Tape Section',
            offset=0.0, offsetType=MIDDLE_SURFACE, offsetField='',
            thicknessAssignment=FROM_SECTION)
    else:
        mirrorPart.SectionAssignment(
            region=mirrorRegion, sectionName='Undulation Section',
            offset=0.0, offsetType=MIDDLE_SURFACE, offsetField='',
            thicknessAssignment=FROM_SECTION)

    mirrorPart.MaterialOrientation(
        region=mirrorRegion, orientationType=SYSTEM, axis=AXIS_3,
        localCsys=None, fieldName='', 
        additionalRotationType=ROTATION_ANGLE, additionalRotationField='',
        angle=objAngle, stackDirection=STACK_3)

    # mesh mirrored part
    mirrorPart.seedPart(size=1.0, deviationFactor=0.1, minSizeFactor=0.1)
    mirrorPart.generateMesh()

    # create the part instance
    instanceName = '{} Instance {}-{}'.format(objType, layer, definePart.counter)
    mirrorInstanceName = '{} Mirror Instance {}-{}'.format(objType, layer, definePart.counter)
    objInstance = laminateAssembly.Instance(
            name=instanceName, part=objPart, dependent=ON)
    laminateAssembly.translate(
            instanceList=(instanceName, ), vector=(0.0, 0.0, layer*cpt))
    mirrorObjInstance = laminateAssembly.Instance(
            name=mirrorInstanceName, part=mirrorPart, dependent=ON)
    laminateAssembly.translate(
            instanceList=(mirrorInstanceName, ), vector=(0.0, 0.0, -layer*cpt))

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
        resinMergedByAngle.append((ang,cascaded_union(resinObjByAngle)))
    for (resinAngle,resinObjAngleSet) in resinMergedByAngle:
        if resinObjAngleSet.geom_type == 'Polygon':
            # angle not important: resin is isotropic
            definePart(resinObjAngleSet, 'Resin',
                       resinAngle, layerNumber) 
        elif resinObjAngleSet.geom_type == 'MultiPolygon':
            for resinPolyObj in resinObjAngleSet:
                definePart(resinPolyObj, 'Resin', resinAngle,
                           layerNumber)

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
            region=passRegion, orientationType=SYSTEM, axis=AXIS_3,
            localCsys=None, fieldName='', additionalRotationType=ROTATION_ANGLE,
            additionalRotationField='', angle=pathAngle, stackDirection=STACK_3)

# IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
# Create cohesive zone interactions between faces of the assembly

# explicit
laminateModel.ContactProperty('Tangential')
laminateModel.interactionProperties['Tangential'].TangentialBehavior(
    formulation=PENALTY, directionality=ISOTROPIC, slipRateDependency=OFF, 
    pressureDependency=OFF, temperatureDependency=OFF, dependencies=0, 
    table=((0.15, ), ), shearStressLimit=None, maximumElasticSlip=FRACTION, 
    fraction=0.005, elasticSlipStiffness=None)
laminateModel.interactionProperties['Tangential'].NormalBehavior(
    pressureOverclosure=HARD, allowSeparation=ON, 
    constraintEnforcementMethod=DEFAULT)
laminateModel.ContactProperty('Cohesive')
laminateModel.interactionProperties['Cohesive'].CohesiveBehavior()
# laminateModel.interactionProperties['Cohesive'].Damage(initTable=((
#     1.0, 2.0, 3.0), ), evolTable=())

# determine contacts
laminateModel.contactDetection(defaultType=CONTACT,
    interactionProperty='Cohesive', nameEachSurfaceFound=OFF,
    createUnionOfMasterSurfaces=ON, createUnionOfSlaveSurfaces=ON)

laminateModel.ContactExp(name='GC', createStepName='Initial')
laminateModel.interactions['GC'].includedPairs.setValuesInStep(
    stepName='Initial', useAllstar=ON)
laminateModel.interactions['GC'].contactPropertyAssignments.appendInStep(
        stepName='Initial', assignments=((GLOBAL, SELF, 'Tangential'), ))
for interact in laminateModel.interactions.values():
    if interact.name == 'GC': 
        continue
    else:
        masterName = interact.master[0]
        slaveName = interact.slave[0]
        masterSurf = laminateAssembly.surfaces[masterName]
        slaveSurf = laminateAssembly.surfaces[slaveName]
        laminateModel.interactions['GC'].contactPropertyAssignments.appendInStep(
            stepName='Initial', assignments=((masterSurf, slaveSurf, 'Cohesive'), ))
        del laminateModel.interactions['{}'.format(interact.name)]


# IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
# Create rigid bodies for four point bending

rollerRadius = 2.5
rollerSketch = laminateModel.ConstrainedSketch(name='Roller', 
    sheetSize=200.0)
rollerSketch.CircleByCenterPerimeter(
        center=(0.0, 0.0), point1=(rollerRadius, 0.0))
rollerPart = laminateModel.Part(name='Roller Part', dimensionality=THREE_D, 
    type=DISCRETE_RIGID_SURFACE)
rollerPart.BaseSolidExtrude(sketch=rollerSketch, depth=50.0)
rollerCells = rollerPart.cells
rollerPart.RemoveCells(cellList = rollerCells[0:1])
rollerRFID = rollerPart.ReferencePoint(point=rollerPart.InterestingPoint(
        edge=rollerPart.edges[0], rule=CENTER)).id
# mesh roller part
rollerPart.seedPart(size=0.71, deviationFactor=0.1, minSizeFactor=0.1)
rollerPart.generateMesh()

# Place rigid bodies in assembly
# roller 1
laminateAssembly.Instance(name='Roller 1', part=rollerPart, dependent=ON)
laminateAssembly.rotate(instanceList=('Roller 1', ), axisPoint=(0.0, 0.0, 0.0), 
    axisDirection=(0.0, 1.0, 0.0), angle=90.0)
laminateAssembly.translate(
        instanceList=('Roller 1', ),
        vector=(-25.0, 25.0, rollerRadius+tLaminate))
# roller 2
laminateAssembly.Instance(name='Roller 2', part=rollerPart, dependent=ON)
laminateAssembly.rotate(instanceList=('Roller 2', ), axisPoint=(0.0, 0.0, 0.0), 
    axisDirection=(0.0, 1.0, 0.0), angle=90.0)
laminateAssembly.translate(
        instanceList=('Roller 2', ),
        vector=(-25.0, -25.0, rollerRadius+tLaminate))
# roller 3
laminateAssembly.Instance(name='Roller 3', part=rollerPart, dependent=ON)
laminateAssembly.rotate(instanceList=('Roller 3', ), axisPoint=(0.0, 0.0, 0.0), 
    axisDirection=(0.0, 1.0, 0.0), angle=90.0)
laminateAssembly.translate(
        instanceList=('Roller 3', ),
        vector=(-25.0, 50.0, -rollerRadius-tLaminate))
# roller 4
laminateAssembly.Instance(name='Roller 4', part=rollerPart, dependent=ON)
laminateAssembly.rotate(instanceList=('Roller 4', ), axisPoint=(0.0, 0.0, 0.0), 
    axisDirection=(0.0, 1.0, 0.0), angle=90.0)
laminateAssembly.translate(
    instanceList=('Roller 4', ),
    vector=(-25.0, -50.0, -rollerRadius-tLaminate))

# IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
# Apply boundary conditions

# identify faces at top and bottom for coupling
tFaces = [inst.faces.getByBoundingBox(-60.0,-76.0,0.34,60.0,76.0,0.36) for inst
                in laminateAssembly.instances.values() 
                if inst.faces.getByBoundingBox(-60.0,-76.0,0.34,60.0,76.0,0.36)]
tFacesSurface = laminateAssembly.Surface(
        side1Faces=tFaces, name='Top Faces')
bFaces = [inst.faces.getByBoundingBox(-60.0,-76.0,-0.36,60.0,76.0,-0.34) for inst
                in laminateAssembly.instances.values() 
                if inst.faces.getByBoundingBox(-60.0,-76.0,-0.36,60.0,76.0,-0.34)]
bFacesSurface = laminateAssembly.Surface(
        side1Faces=bFaces, name='Bottom Faces')

# explicit
laminateModel.SmoothStepAmplitude(name='Smoothing Amplitude',
    timeSpan=STEP, data=((0.0, 0.0), (1e-05, 1.0)))

rfPoint1 = laminateAssembly.instances['Roller 1'].referencePoints[3]
rfPoint2 = laminateAssembly.instances['Roller 2'].referencePoints[3]
rfPoint3 = laminateAssembly.instances['Roller 3'].referencePoints[3]
rfPoint4 = laminateAssembly.instances['Roller 4'].referencePoints[3]

topRollerRegion = laminateAssembly.Set(referencePoints=(rfPoint1, rfPoint2, ),
        name='Top Rollers')
laminateModel.DisplacementBC(name='Top Roller BC', createStepName='Initial', 
    region=topRollerRegion, u1=SET, u2=SET, u3=UNSET, ur1=SET, ur2=SET, ur3=SET, 
    amplitude=UNSET, distributionType=UNIFORM, fieldName='', 
    localCsys=None)
bottomRollerRegion = laminateAssembly.Set(referencePoints=(rfPoint3, rfPoint4, ),
        name='Bottom Rollers')
laminateModel.EncastreBC(name='Bottom Roller BC', createStepName='Initial', 
    region=bottomRollerRegion, localCsys=None)

laminateModel.DisplacementBC(name='Top Roller Displacement', 
    createStepName='Loading Step', region=topRollerRegion, u1=0.0, u2=0.0, 
    u3=-0.025, ur1=0.0, ur2=0.0, ur3=0.0, amplitude='Smoothing Amplitude', 
    fixed=OFF, distributionType=UNIFORM, fieldName='', localCsys=None)


