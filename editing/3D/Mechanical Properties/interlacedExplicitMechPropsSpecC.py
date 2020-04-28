from sys import path
githubPath = r"C:\Users\rutge\Documents\GitHub"
path.append(r'C:\Python27\Lib\site-packages')
path.append(githubPath + '\\interlaced_model_creation\\editing')
from abaqus import *
from abaqusConstants import *
import regionToolset
import sigc as sigc
import tapePlacement as tp
from itertools import combinations
import analyticStiffness
import math
import ast
from shapely.geometry import Polygon

session.viewports['Viewport: 1'].setValues(displayedObject=None)

# IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
# Define model parameters

specimenType = 'C'
laminateAngles = (0, 90)  # define angles of tapes
tapeWidth = 15.0
tw = (tapeWidth, ) * len(laminateAngles)  # tape widths
# number of gaps between tapes in interlacing pattern
tapeSpace = 2
# cured ply thickness e.g. 0.18125
cpt = 0.18
# ratio of undulation amplitude to length e.g. 0.090625
undulationRatio = 0.18
uw = cpt / undulationRatio

# define specimen dimensions
if specimenType == 'B' or 'C':
    xMin = -45.0 / 2.0
    xMax = -xMin
else:
    xMin = -25.0
    xMax = 25.0
yMin = -75.0
yMax = -yMin
specimenHeight = yMax - yMin
specimenWidth = xMax - xMin
numLayers = len(laminateAngles) * 2.0  # symmetric
laminateThickness = numLayers * cpt
zMax = laminateThickness / 2.0
zMin = -zMax
RVEPolygon = Polygon([(xMin, yMin), (xMax, yMin), (xMax, yMax), (xMin, yMax)])
dimensions = [xMin, yMin, zMin, xMax, yMax, zMax]
meshSize = 1.0

# create grid using shapely
partGrid = sigc.createGrids(tapeAngles=laminateAngles, tapeWidths=tw,
                            undulationWidth=uw, sample=RVEPolygon)
# identify parts in grid
tapePaths = tp.laminateCreation(
    grid=partGrid, tapeAngles=laminateAngles, tapeWidths=tw,
    tapeSpacing=tapeSpace, undulationWidth=uw)

# IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
# Define model name

ratioString = str(undulationRatio).replace('.', '_')
modelName = '{}-{}'.format(specimenType, ratioString)
mdb.Model(name=modelName, modelType=STANDARD_EXPLICIT)
laminateModel = mdb.models[modelName]
laminateAssembly = laminateModel.rootAssembly

# set the work directory
wd = 'C:\\Workspace\\3D_MechProps\\{}'.format(modelName)
if not os.path.exists(wd):
    os.makedirs(wd)
os.chdir(wd)

# IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
# Material input data
E11 = 146.8
E22 = E33 = 11.6
nu12 = nu13 = 0.289
nu23 = 0.298
G12 = G13 = 6.47
G23 = 4.38
Xt = 2.61
Xc = 1.759
Yt = 0.055
Yc = 0.285
Sl = 0.105
alpha0 = 53.0
G1Plus = 0.1
G1Minus = 0.1
G2Plus = 0.00075
G2Minus = 0.0025
G6 = 0.0035
density = 1.59e-06

# Create Abaqus materials

tapeMaterial = laminateModel.Material(name='Tape')
tapeMaterial.Density(table=((density, ), ))
tapeMaterial.Depvar(deleteVar=20, n=24)
tapeMaterial.UserMaterial(
    mechanicalConstants=(E11, E22, E33, nu12, nu13, nu23, G12, G13,
                         G23, Xt, Xc, Yt, Yc, Sl, alpha0, G1Plus,
                         G1Minus, G2Plus, G2Minus, G6))
tapeSection = laminateModel.HomogeneousSolidSection(name='Tape Section',
                                                    material='Tape')

resinMaterial = laminateModel.Material(name='Resin')
resinMaterial.Elastic(type=ENGINEERING_CONSTANTS,
                      table=((7.47, 7.47, 7.47,
                              0.32, 0.32, 0.32,
                              3.5, 3.5, 3.5), ))
resinMaterial.Density(table=((1.1e-06, ), ))
resinSection = laminateModel.HomogeneousSolidSection(name='Resin Section',
                                                     material='Resin')

# Determine stiffness of undulation regions using analytical model
# Adapted from Chou et al. 1972

interfaceAngleCombos = [(0, ang) for ang in laminateAngles if ang != 0]
n = 200  # for analytic model discretization
for combo in interfaceAngleCombos:
    matName = 'Undulation {}'.format(combo)
    O1 = math.radians(combo[1])  # in-plane angle of bottom ply
    O2 = math.radians(combo[0])  # in-plane angle of undulating ply
    Clamina = analyticStiffness.CFromConstants(E11, E22, nu12, nu23, G12, G23)
    Vfrac, a = analyticStiffness.undulationGeometry(cpt, uw, n)
    CIsostrain, CIsostress = analyticStiffness.determineStiffness(
        Clamina, a, O1, O2, Vfrac, n)
    engConstants = analyticStiffness.engineeringConstants(CIsostrain)
    matProps = engConstants + [2.61, 1.759, 0.055, 0.285, 0.105, 53.0, 0.1,
                               0.1, 0.00075, 0.0025, 0.0035]
    undMaterial = laminateModel.Material(name=matName)
    undMaterial.Depvar(deleteVar=20, n=24)
    undMaterial.UserMaterial(mechanicalConstants=matProps)
    undMaterial.Density(table=((density, ), ))
    undulationSection = laminateModel.HomogeneousSolidSection(
        name=matName + ' Section', material=matName)

# IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
# Create step

# create step
laminateModel.ExplicitDynamicsStep(name='Loading Step', previous='Initial',
                                   timePeriod=20.0,
                                   massScaling=((SEMI_AUTOMATIC, MODEL,
                                                 THROUGHOUT_STEP, 0.0, 1e-06,
                                                 BELOW_MIN, 1, 0, 0.0, 0.0, 0,
                                                 None), ))

# Modify field output request
fieldOutput = laminateModel.fieldOutputRequests['F-Output-1']
fieldOutput.setValues(variables=('S', 'LE', 'U', 'V', 'A', 'RF', 'CSTRESS',
                                 'CSDMG', 'SDV', 'STATUS'), numIntervals=50)


# IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
# Generate parts

startNodeNum = 1
startElemNum = 1


# This function is used to create parts using the geometric info from Shapely
def definePart(obj, objectID):
    uniqueID = '{}-{}'.format(obj.layer, objectID)
    x0, y0 = obj.exterior.xy
    abqCoords = zip(x0, y0)

    profileSketch = laminateModel.ConstrainedSketch(
        name='Sketch {}'.format(uniqueID), sheetSize=200.0)
    for ind in range(len(abqCoords) - 1):
        profileSketch.Line(point1=(abqCoords[ind][0], abqCoords[ind][1]),
                           point2=(abqCoords[ind + 1][0],
                                   abqCoords[ind + 1][1]))

    # extrude profile to create part
    objPart = laminateModel.Part(name='Part {}'.format(uniqueID),
                                 dimensionality=THREE_D, type=DEFORMABLE_BODY)
    objPart.BaseSolidExtrude(sketch=profileSketch, depth=cpt)

    objCells = objPart.cells  # all cells within part (only 1 cell)
    objRegion = regionToolset.Region(cells=objCells[:])

    if obj.objectType == 'Tape':
        objPart.SectionAssignment(
            region=objRegion, sectionName='Tape Section',
            offset=0.0, offsetType=MIDDLE_SURFACE, offsetField='',
            thicknessAssignment=FROM_SECTION)
    elif obj.objectType == 'Undulation':
        interfaceAngle = abs(obj.angle[0] - obj.angle[1])
        secAngle = (0, interfaceAngle)
        secName = 'Undulation {} Section'.format(secAngle)
        objPart.SectionAssignment(
            region=objRegion, sectionName=secName,
            offset=0.0, offsetType=MIDDLE_SURFACE, offsetField='',
            thicknessAssignment=FROM_SECTION)
    else:
        objPart.SectionAssignment(
            region=objRegion, sectionName='Resin Section',
            offset=0.0, offsetType=MIDDLE_SURFACE, offsetField='',
            thicknessAssignment=FROM_SECTION)

    objPart.MaterialOrientation(
        region=objRegion, orientationType=SYSTEM, axis=AXIS_3,
        localCsys=None, fieldName='', additionalRotationType=ROTATION_ANGLE,
        additionalRotationField='', angle=obj.angle[-1],
        stackDirection=STACK_3)

    # create the part instance
    instanceName = 'Instance {}'.format(uniqueID)
    objInstance = laminateAssembly.Instance(name=instanceName, part=objPart,
                                            dependent=ON)
    laminateAssembly.translate(instanceList=(instanceName, ),
                               vector=(0.0, 0.0, (obj.layer - 1) * cpt))

    # mesh the part instance
    global startElemNum, startNodeNum
    objPart.setValues(startNodeLabel=startNodeNum, startElemLabel=startElemNum)
    objPart.seedPart(size=meshSize)
    objPart.generateMesh()
    startElemNum += len(objPart.elements)
    startNodeNum += len(objPart.nodes)


# ----------------------------------------------------------------------------
# In the next section the Shapely geometry objects are iterated over and
# definePart() is called to create the parts

for layerNumber in range(1, len(partGrid) + 1):
    for objID, obj in partGrid[layerNumber].iteritems():
        definePart(obj, objID)

# this section merges all the parts in a tapepath into one part
for (pathNumber, tapePath) in enumerate(tapePaths):
    if len(tapePath) > 1:
        layer, gridID = tapePath[0].split('-')
        pathAngle = partGrid[int(layer)][int(gridID)].angle[-1]
        instanceList = [laminateAssembly.instances['Instance ' + str(inst)]
                        for inst in tapePath]
        instanceName = 'Pass {}'.format(pathNumber)
        laminateAssembly.InstanceFromBooleanMerge(name=instanceName,
                                                  instances=instanceList,
                                                  keepIntersections=ON,
                                                  originalInstances=DELETE,
                                                  mergeNodes=ALL,
                                                  nodeMergingTolerance=1e-06,
                                                  domain=BOTH)
        passPart = laminateModel.parts[instanceName]
        passRegion = regionToolset.Region(cells=passPart.cells[:])
        passPart.MaterialOrientation(region=passRegion, orientationType=SYSTEM,
                                     axis=AXIS_3, localCsys=None, fieldName='',
                                     additionalRotationType=ROTATION_ANGLE,
                                     additionalRotationField='',
                                     angle=pathAngle, stackDirection=STACK_3)
        passPart.setValues(startNodeLabel=startNodeNum,
                           startElemLabel=startElemNum)
        passPart.seedPart(size=meshSize)
        passPart.generateMesh()
        startElemNum += len(passPart.elements)
        startNodeNum += len(passPart.nodes)

allInstanceKeys = laminateAssembly.instances.keys()
laminateAssembly.LinearInstancePattern(instanceList=allInstanceKeys, 
                                       direction1=(0.0, 0.0, 1.0),
                                       direction2=(0.0, 1.0, 0.0), number1=2, 
                                       number2=1, spacing1=0.36, spacing2=150.0)

# IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
# Create cohesive zone interactions between faces of the assembly
# Calculate CZM parameters according to Turon 2007 (sort of...)

alpha = 50
E_33 = tapeMaterial.userMaterial.mechanicalConstants[2]
K = (alpha * E_33) / cpt

# define contact properties
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
laminateModel.interactionProperties['Cohesive'].CohesiveBehavior(
    defaultPenalties=OFF, table=((K, K, K), ))
laminateModel.interactionProperties['Cohesive'].Damage(
    criterion=QUAD_TRACTION, initTable=((0.055, 0.09, 0.09), ),
    useEvolution=ON, evolutionType=ENERGY, useMixedMode=ON,
    mixedModeType=POWER_LAW, exponent=1.0,
    evolTable=((3.52e-4, 3.52e-4, 3.52e-4), ))

# determine contacts
laminateModel.contactDetection(defaultType=CONTACT,
    interactionProperty='Cohesive', nameEachSurfaceFound=OFF,
    createUnionOfMasterSurfaces=ON, createUnionOfSlaveSurfaces=ON,
    separationTolerance=cpt/200.0)

# convert individual contact interactions to a single interaction with
# each contact defined as an 'Individual Property Assignment'
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
            stepName='Initial',
            assignments=((masterSurf, slaveSurf, 'Cohesive'), ))
        del laminateModel.interactions['{}'.format(interact.name)]

# IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
# Apply boundary conditions

tol = 0.001
# identify faces at top and bottom for coupling
tFaces = [inst.faces.getByBoundingBox(xMin - tol, yMax - tol, zMin - tol,
                                      xMax + tol, yMax + tol, zMax + tol)
          for inst in laminateAssembly.instances.values()
          if inst.faces.getByBoundingBox(xMin - tol, yMax - tol, zMin - tol,
                                         xMax + tol, yMax + tol, zMax + tol)]
tFacesSurface = laminateAssembly.Surface(side1Faces=tFaces, name='Top Faces')
bFaces = [inst.faces.getByBoundingBox(xMin - tol, yMin - tol, zMin - tol,
                                      xMax + tol, yMin + tol, zMax + tol)
          for inst in laminateAssembly.instances.values()
          if inst.faces.getByBoundingBox(xMin - tol, yMin - tol, zMin - tol,
                                         xMax + tol, yMin + tol, zMax + tol)]
bFacesSurface = laminateAssembly.Surface(side1Faces=bFaces, name='Bottom Faces')
symFaces = [inst.faces.getByBoundingBox(xMin - tol, yMin - tol, -tol,
                                        xMax + tol, yMax + tol, tol)
            for inst in laminateAssembly.instances.values()
            if inst.faces.getByBoundingBox(xMin - tol, yMin - tol, -tol,
                                           xMax + tol, yMax + tol, tol)]
symFacesSet = laminateAssembly.Set(faces=symFaces, name='Symmetry Faces')

# create reference points for coupling
rfPoint1Id = laminateAssembly.ReferencePoint(point=(0.0, yMax + 5.0, 0.0)).id
rfPoint1 = laminateAssembly.referencePoints[rfPoint1Id]
rfPoint1Region = laminateAssembly.Set(referencePoints=(rfPoint1,),
                                      name='Coupling Reference Point 1')
rfPoint2Id = laminateAssembly.ReferencePoint(point=(0.0, yMin - 5.0, 0.0)).id
rfPoint2 = laminateAssembly.referencePoints[rfPoint2Id]
rfPoint2Region = laminateAssembly.Set(referencePoints=(rfPoint2,),
                                      name='Coupling Reference Point 2')

laminateModel.Coupling(name='Top Coupling', controlPoint=rfPoint1Region,
                       surface=tFacesSurface, influenceRadius=WHOLE_SURFACE,
                       couplingType=KINEMATIC, localCsys=None, u1=ON, u2=ON,
                       u3=ON, ur1=ON, ur2=ON, ur3=ON)

laminateModel.Coupling(name='Bottom Coupling', controlPoint=rfPoint2Region,
                       surface=bFacesSurface, influenceRadius=WHOLE_SURFACE,
                       couplingType=KINEMATIC, localCsys=None, u1=ON, u2=ON,
                       u3=ON, ur1=ON, ur2=ON, ur3=ON)

# explicit
laminateModel.SmoothStepAmplitude(name='Smoothing Amplitude', timeSpan=STEP,
                                  data=((0.0, 0.0), (1e-05, 1.0)))
laminateModel.DisplacementBC(name='Bottom Surface BC',
                             createStepName='Loading Step',
                             region=rfPoint2Region, u1=0.0, u2=0.0, u3=0.0,
                             ur1=0.0, ur2=0.0, ur3=0.0, amplitude=UNSET,
                             fixed=OFF, distributionType=UNIFORM, fieldName='',
                             localCsys=None)
laminateModel.VelocityBC(name='Top Surface BC', createStepName='Loading Step',
                         region=rfPoint1Region, v1=UNSET, v2=0.05, v3=UNSET,
                         vr1=UNSET, vr2=UNSET, vr3=UNSET,
                         amplitude='Smoothing Amplitude', localCsys=None,
                         distributionType=UNIFORM, fieldName='')
laminateModel.ZsymmBC(name='Symmetry', createStepName='Loading Step',
                      region=symFacesSet, localCsys=None)
