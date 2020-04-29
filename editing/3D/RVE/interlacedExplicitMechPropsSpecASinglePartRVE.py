from sys import path
path.append(r'C:\Python27\Lib\site-packages')
path.append(r"\\arran.sms.ed.ac.uk\home\s1342398\GitHub\interlaced_model_creation\editing")
from shapely.geometry import Point, Polygon, LineString
from shapely import affinity
from abaqus import *
from abaqusConstants import *
from visualization import *
import regionToolset
import sigc as sigc
import tapePlacement as tp
import analyticStiffness
import math
import os
import mesh
import ast

session.viewports['Viewport: 1'].setValues(displayedObject=None)

# IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
# Define model parameters
userInput = getInputs(fields=(('Laminate Angles (deg):', '(0, 90)'),
                              ('Tape Width (mm):', '25'),
                              ('Tape Spacing (int):', '1'),
                              ('Cured Ply Thickness (mm):', '0.18'),
                              ('Undulation Ratio (float):', '0.09')),
                      label='Please provide the following information',
                      dialogTitle='Model Parameters')

laminateAngles = ast.literal_eval(userInput[0])  # define angles of tapes
tapeWidth = float(userInput[1])
tw = (tapeWidth, ) * len(laminateAngles)  # tape widths
# number of gaps between tapes in interlacing pattern
tapeSpace = int(userInput[2])
# cured ply thickness e.g. 0.18125
cpt = float(userInput[3])
# ratio of undulation amplitude to length e.g. 0.090625
undulationRatio = float(userInput[4])
uw = cpt / undulationRatio

# define RVE dimensions
xMin = yMin = -(tw[0] / 2.0)
xMax = yMax = xMin + (tapeSpace + 1) * (tw[0])
specimenHeight = yMax - yMin
specimenWidth = xMax - xMin
xMid = yMid = xMin + (specimenWidth / 2.0)
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

ratioString = userInput[4].replace('.', '_')
angleString = '-'.join([str(x) for x in laminateAngles])
modelName = '{}--{}-{}-{}'.format(angleString, int(tw[0]), int(tapeSpace),
                                  ratioString)
mdb.Model(name=modelName, modelType=STANDARD_EXPLICIT)
laminateModel = mdb.models[modelName]
laminateAssembly = laminateModel.rootAssembly

# set the work directory
wd = 'C:\\Workspace\\3D_RVE\\{}'.format(modelName)
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
                                   timePeriod=60.0,
                                   massScaling=((SEMI_AUTOMATIC, MODEL,
                                                 THROUGHOUT_STEP, 0.0, 1e-06,
                                                 BELOW_MIN, 1, 0, 0.0, 0.0, 0,
                                                 None), ))
# Modify field output request
fieldOutput = laminateModel.fieldOutputRequests['F-Output-1']
fieldOutput.setValues(variables=('S', 'LE', 'U', 'V', 'A', 'RF', 'SDV',
                                 'STATUS'), numIntervals=50)

# IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
# Generate parts

# Create sketch
specimenSketch = laminateModel.ConstrainedSketch(name='Specimen Sketch',
                                                 sheetSize=200.0)
specimenSketch.rectangle(point1=(xMin, yMin), point2=(xMax, yMax))
specimenPart = laminateModel.Part(name='Specimen',
                                  dimensionality=THREE_D,
                                  type=DEFORMABLE_BODY)
specimenPart.BaseSolidExtrude(sketch=specimenSketch, depth=zMax)

# Through thickness partition into plies
plyDatumPlaneIds = []
for k in range(1, len(laminateAngles)):
    plyDatumPlaneIds.append(specimenPart.DatumPlaneByPrincipalPlane(
        principalPlane=XYPLANE, offset=cpt * k).id)
    specimenCells = specimenPart.cells
    specimenPart.PartitionCellByDatumPlane(
        datumPlane=specimenPart.datums[plyDatumPlaneIds[k - 1]],
        cells=specimenCells)


def partitionSpecimen(tapeAngles, tapeWidths, undulationWidth=1.0):
    '''
    Function to partition the specimen into a grid of polygons depending on the
    tape angles of the laminate.

    '''
    # define mirror point (used to mirror the Polygon boundaries)
    mirrorPoint = Point([(0.0, 0.0), (0.0, 0.0)])

    tapes = zip(tapeAngles, tapeWidths)
    partitionLines = []
    uw = undulationWidth
    for tape in tapes:
        a = tape[0]  # angle
        w = tape[1]  # width
        if a == 90:
            offset = w
            maxOffset = xMax - w / 2.0
            numberOffsets = int(maxOffset / offset)
            offsetList = [offset * k for k
                          in range(0, numberOffsets + 1)]
            for ofs in offsetList:
                tapeLineCoords = [((w / 2.0) - uw + ofs, 100.0),
                                  ((w / 2.0) - uw + ofs, -100.0)]
                resinLineCoords1 = ([((w / 2.0) + ofs, 100.0),
                                     ((w / 2.0) + ofs, -100.0)])
                resinLineCoords2 = ([((w / 2.0) + ofs + uw, 100.0),
                                     ((w / 2.0) + ofs + uw, -100.0)])
                tapeLine = LineString(tapeLineCoords)
                resinLine1 = LineString(resinLineCoords1)
                resinLine2 = LineString(resinLineCoords2)
                partitionLines.extend([tapeLine, resinLine1, resinLine2])
            reflectedLines = [affinity.scale(line, xfact=-1,
                                             origin=mirrorPoint)
                              for line in partitionLines]
            partitionLines = partitionLines + reflectedLines

        else:
            offset = w / math.cos(math.radians(a))
            maxOffset = (yMax - (w / 2.0) * math.cos(math.radians(a))
                         + 50.0 * math.tan(math.radians(a)))
            numberOffsets = int(maxOffset / offset)
            offsetList = [offset * k for k in range(0, numberOffsets + 1)]
            for ofs in offsetList:
                tapeLineCoords = [(-100.0, (w / 2.0) - uw + ofs),
                                  (100.0, (w / 2.0) - uw + ofs)]
                resinLineCoords1 = ([(-100.0, (w / 2.0) + ofs),
                                     (100.0, (w / 2.0) + ofs)])
                resinLineCoords2 = ([(-100.0, (w / 2.0) + ofs + uw),
                                     (100.0, (w / 2.0) + ofs + uw)])
                rotPoint = Point([(0.0, ofs), (0.0, ofs)])
                tapeLine = affinity.rotate(
                    LineString(tapeLineCoords), a, rotPoint)
                resinLine1 = affinity.rotate(
                    LineString(resinLineCoords1), a, rotPoint)
                resinLine2 = affinity.rotate(
                    LineString(resinLineCoords2), a, rotPoint)
                partitionLines.extend([tapeLine, resinLine1, resinLine2])
            reflectedLines = [affinity.rotate(line, 180.0, origin=mirrorPoint)
                              for line in partitionLines]
            partitionLines = partitionLines + reflectedLines

    for line in partitionLines:
        lineCoords = [list(tup) + [0.0] for tup in list(line.coords)]
        finalCoord = [[lineCoords[-1][0], lineCoords[-1][1], 1.0]]
        lineCoords3D = lineCoords + finalCoord
        dPointIDs = []
        for coordinates in lineCoords3D:
            dPointIDs.append(specimenPart.DatumPointByCoordinate(
                coords=coordinates).id)
        dPoint1 = specimenPart.datums[dPointIDs[0]]
        dPoint2 = specimenPart.datums[dPointIDs[1]]
        dPoint3 = specimenPart.datums[dPointIDs[2]]
        dpbtp = specimenPart.DatumPlaneByThreePoints(point1=dPoint1,
                                                     point2=dPoint2,
                                                     point3=dPoint3).id
        specimenCells = specimenPart.cells
        try:
            specimenPart.PartitionCellByDatumPlane(
                datumPlane=specimenPart.datums[dpbtp], cells=specimenCells)
        except:
            continue


partitionSpecimen(tapeAngles=laminateAngles, tapeWidths=tw, undulationWidth=uw)

# IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
# Mesh

specimenPart.seedPart(size=meshSize)
specimenPart.generateMesh()
elemType1 = mesh.ElemType(elemCode=C3D8R, elemLibrary=EXPLICIT, 
                          kinematicSplit=AVERAGE_STRAIN,
                          secondOrderAccuracy=OFF, hourglassControl=ENHANCED,
                          distortionControl=DEFAULT)
elemType2 = mesh.ElemType(elemCode=C3D6, elemLibrary=EXPLICIT)
elemType3 = mesh.ElemType(elemCode=C3D4, elemLibrary=EXPLICIT)
allCellsSet = specimenPart.Set(cells=specimenPart.cells, name='All Cells')
specimenPart.setElementType(regions=allCellsSet,
                            elemTypes=(elemType1, elemType2, elemType3))

# IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
# Assign properties

for layer, polyDict in partGrid.iteritems():
    for objID, obj in polyDict.iteritems():
        objCentroid = obj.centroid.coords
        zCoord = (2 * layer - 1) / 2.0 * cpt
        selectedCells = specimenCells.findAt(((objCentroid[0][0],
                                               objCentroid[0][1], zCoord), ),
                                             ((objCentroid[0][0],
                                               objCentroid[0][1], -zCoord), ))
        objRegion = regionToolset.Region(cells=selectedCells)
        objType = obj.objectType
        objAngle = obj.angle[-1]
        if objType == 'Undulation':
            interfaceAngle = abs(obj.angle[0] - obj.angle[1])
            secAngle = (0, interfaceAngle)
            secName = 'Undulation {} Section'.format(secAngle)
        else:
            secName = '{} Section'.format(objType)
        specimenPart.SectionAssignment(
            region=objRegion, sectionName=secName,
            offset=0.0, offsetType=MIDDLE_SURFACE, offsetField='',
            thicknessAssignment=FROM_SECTION)
        specimenPart.MaterialOrientation(
            region=objRegion, orientationType=SYSTEM, axis=AXIS_3,
            localCsys=None, fieldName='',
            additionalRotationType=ROTATION_ANGLE, additionalRotationField='',
            angle=objAngle, stackDirection=STACK_3)

# create the part instance
specimenInstance = laminateAssembly.Instance(name='Specimen',
                                             part=specimenPart, dependent=ON)

# IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
# Apply boundary conditions

tol = 0.001
# identify faces at top and bottom for coupling
tFaces = specimenInstance.faces.getByBoundingBox(xMin - tol, yMax - tol,
                                                 zMin - tol, xMax + tol,
                                                 yMax + tol, zMax + tol)
tFacesSurface = laminateAssembly.Surface(side1Faces=tFaces, name='Top Faces')
bFaces = specimenInstance.faces.getByBoundingBox(xMin - tol, yMin - tol,
                                                 zMin - tol, xMax + tol,
                                                 yMin + tol, zMax + tol)
bFacesSurface = laminateAssembly.Surface(side1Faces=bFaces,
                                         name='Bottom Faces')
symFaces = specimenInstance.faces.getByBoundingBox(xMin - tol, yMin - tol,
                                                   -tol, xMax + tol,
                                                   yMax + tol, tol)
symFacesSet = laminateAssembly.Set(faces=symFaces, name='z Symmetry Faces')
rFaces = specimenInstance.faces.getByBoundingBox(xMax - tol, yMin - tol,
                                                 zMin - tol, xMax + tol,
                                                 yMax + tol, zMax + tol)
rFacesSet = laminateAssembly.Set(faces=rFaces, name='x Symmetry Faces')

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
                       couplingType=KINEMATIC, localCsys=None, u1=OFF, u2=ON,
                       u3=OFF, ur1=ON, ur2=ON, ur3=ON)

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
                         vr1=0.0, vr2=0.0, vr3=0.0,
                         amplitude='Smoothing Amplitude', localCsys=None,
                         distributionType=UNIFORM, fieldName='')
laminateModel.ZsymmBC(name='ZSymmetry', createStepName='Loading Step',
                      region=symFacesSet, localCsys=None)
laminateModel.XsymmBC(name='XSymmetry', createStepName='Loading Step',
                      region=rFacesSet, localCsys=None)

# IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
# Apply boundary conditions

# Create job and write input

# Define path to the material subroutine file (use a raw string)
# vumatPath = r'/home/s1342398/abaqus/catalanotti_FC.f'

# mdb.Job(name=modelName, model=modelName, description='', type=ANALYSIS,
#         atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90,
#         memoryUnits=PERCENTAGE, explicitPrecision=DOUBLE_PLUS_PACK,
#         nodalOutputPrecision=FULL, echoPrint=OFF, modelPrint=OFF,
#         contactPrint=OFF, historyPrint=OFF, userSubroutine=vumatPath,
#         scratch='', resultsFormat=ODB, parallelizationMethodExplicit=DOMAIN,
#         numDomains=4, activateLoadBalancing=False,
#         multiprocessingMode=DEFAULT, numCpus=4)
# mdb.jobs[modelName].submit(consistencyChecking=OFF)

# Save model
# mdb.saveAs(pathName=modelName)

# IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
# Post processing

# Open the Output Database for the current Job
# odb = openOdb(path='{}.odb'.format(modelName))

# frame = odb.steps['Loading Step'].frames[-1]
# rfFieldOutput = frame.fieldOutputs['RF']
# uFieldOutput = frame.fieldOutputs['U']
# masterNode1 = odb.rootAssembly.nodeSets['MASTERNODE1']
# masterNode2 = odb.rootAssembly.nodeSets['MASTERNODE2']
# rfMasterNode1 = rfFieldOutput.getSubset(region=masterNode1)
# uMasterNode1 = uFieldOutput.getSubset(region=masterNode1)

# rf1 = rfMasterNode1.values[0].dataDouble[0]
# u1 = uMasterNode1.values[0].dataDouble[0]

# E11 = ((rf1/(laminateThickness*specimenHeight)) / (u1/specimenWidth))
# print E11
