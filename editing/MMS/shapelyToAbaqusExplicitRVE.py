from sys import path
path.append(r'C:\Python27\Lib\site-packages')
path.append(r"C:\Users\rutge\Documents\GitHub\interlaced_model_creation\editing\MMS")
from shapely.geometry import Point, Polygon, LineString
from shapely import affinity
from abaqus import *
from abaqusConstants import *
import displayGroupMdbToolset as dgm
import regionToolset
import sigc as sigc
import tapePlacement as tp
import analyticStiffness
import math
from periodicBC_2D import periodicBC
import mesh


session.viewports['Viewport: 1'].setValues(displayedObject=None)
'''
Changelog:
- modified script to create shell model with undulations rather than solid
  element model.
- included function to create datum planes and partition the sample part
  into region for material assignment at the Gauss point level.
'''
# IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
# Define model parameters

cpt = 0.175  # cured ply thickness
laminateAngles = (0, 90)  # define angles of tapes in laminate
tw = (25, 25)  # tape widths
uw = 1.0
n = 200  # for analytic model discretization
numLayers = len(laminateAngles)
tapeSpace = 1

xMin = yMin = -(tw[0]/2.0)
xMax = yMax = xMin + (tapeSpace+1)*(tw[0])

RVEPolygon = Polygon(
        [(xMin, yMin), (xMax, yMin), (xMax, yMax), (xMin, yMax)])

partTypes = ['Tape', 'Resin', 'Undulation']
# create grid using shapely
partGrid = sigc.createGrids(tapeAngles=laminateAngles, tapeWidths=tw,
                            undulationWidth=uw, sample=RVEPolygon)
# identify parts in grid
tapePaths = tp.laminateCreation(
    grid=partGrid, tapeAngles=laminateAngles, tapeWidths=tw,
    tapeSpacing=tapeSpace, undulationWidth=uw)

# IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
# Define model name

modelName = 'MMS_Interlaced_RVE'
mdb.Model(name=modelName, modelType=STANDARD_EXPLICIT)
laminateModel = mdb.models[modelName]
laminateAssembly = laminateModel.rootAssembly

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
tapeSection = laminateModel.HomogeneousSolidSection(
        name='Tape Section', material='Tape')

resinMaterial = laminateModel.Material(name='Resin')
resinMaterial.Elastic(type=ENGINEERING_CONSTANTS,
                      table=((7.47, 7.47, 7.47,
                              0.32, 0.32, 0.32,
                              3.5, 3.5, 3.5), ))
resinMaterial.Density(table=((1.1e-06, ), ))
resinSection = laminateModel.HomogeneousSolidSection(
        name='Resin Section', material='Resin')

# Determine stiffness of undulation regions using analytical model
# Adapted from Chou et al. 1972

interfaceAngleCombos = [(0, ang) for ang in laminateAngles if ang != 0]
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
# Step
laminateModel.StaticLinearPerturbationStep(name='Loading Step',
                                           previous='Initial')

laminateModel.fieldOutputRequests['F-Output-1'].setValues(variables=(
    'S', 'LE', 'U', 'RF', 'CF', 'IVOL', 'STH'))

# IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
# Generate parts

# Create sketch
specimenSketch = laminateModel.ConstrainedSketch(name='Specimen Sketch', 
    sheetSize=200.0)
specimenSketch.rectangle(point1=(xMin, yMin), point2=(xMax, yMax))
specimenPart = laminateModel.Part(name='Specimen',
                                  dimensionality=THREE_D,
                                  type=DEFORMABLE_BODY)
specimenPart.BaseShell(sketch=specimenSketch)
specimenFaces = specimenPart.faces

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
            maxOffset = 50.0-w/2.0
            numberOffsets = int(maxOffset/offset)
            offsetList = [offset*k for k
                          in range(0, numberOffsets+1)]
            for os in offsetList:
                tapeLineCoords = [((w/2.0)-uw+os, 100.0),
                                  ((w/2.0)-uw+os, -100.0)]
                resinLineCoords1 = ([((w/2.0)+os, 100.0),
                                     ((w/2.0)+os, -100.0)])
                resinLineCoords2 = ([((w/2.0)+os+uw, 100.0),
                                     ((w/2.0)+os+uw, -100.0)])
                tapeLine = LineString(tapeLineCoords)
                resinLine1 = LineString(resinLineCoords1)
                resinLine2 = LineString(resinLineCoords2)
                partitionLines.extend([tapeLine, resinLine1, resinLine2])
            reflectedLines = [
                    affinity.scale(line, xfact=-1, origin=mirrorPoint)
                    for line in partitionLines]
            partitionLines = partitionLines + reflectedLines

        else:
            offset = w/math.cos(math.radians(a))
            maxOffset = (75.0 - (w/2.0)*math.cos(math.radians(a))
                         + 50.0*math.tan(math.radians(a)))
            numberOffsets = int(maxOffset/offset)
            offsetList = [offset*k for k in range(0, numberOffsets+20)]
            for os in offsetList:
                tapeLineCoords = [(-100.0, (w/2.0)-uw+os),
                                  (100.0, (w/2.0)-uw+os)]
                resinLineCoords1 = (
                        [(-100.0, (w/2.0)+os), (100.0, (w/2.0)+os)])
                resinLineCoords2 = (
                        [(-100.0, (w/2.0)+os+uw), (100.0, (w/2.0)+os+uw)])
                rotPoint = Point([(0.0, os), (0.0, os)])
                tapeLine = affinity.rotate(
                    LineString(tapeLineCoords), a, rotPoint)
                resinLine1 = affinity.rotate(
                    LineString(resinLineCoords1), a, rotPoint)
                resinLine2 = affinity.rotate(
                    LineString(resinLineCoords2), a, rotPoint)
                partitionLines.extend([tapeLine, resinLine1, resinLine2])
            reflectedLines = [
                    affinity.rotate(line, 180.0, origin=mirrorPoint)
                    for line in partitionLines]
            partitionLines = partitionLines + reflectedLines

    for line in partitionLines:
        lineCoords = [list(tup)+[0.0] for tup in list(line.coords)]
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
        specimenFaces = specimenPart.faces
        try:
            specimenPart.PartitionFaceByDatumPlane(
                datumPlane=specimenPart.datums[dpbtp], faces=specimenFaces)
        except:
            continue

partitionSpecimen(tapeAngles=laminateAngles, tapeWidths=tw, undulationWidth=uw)

for objID, obj in partGrid[1].iteritems():
    objCentroid = obj.centroid.coords
    selectedFace = specimenFaces.findAt(((objCentroid[0][0],
                                          objCentroid[0][1], 0.0), ))
    region = regionToolset.Region(faces=selectedFace)
    compositeLayup = specimenPart.CompositeLayup(
        name='CompositeLayup-{}'.format(objID), description='',
        elementType=SHELL, offsetType=MIDDLE_SURFACE, symmetric=True, 
        thicknessAssignment=FROM_SECTION)
    compositeLayup.Section(preIntegrate=OFF, integrationRule=SIMPSON, 
        thicknessType=UNIFORM, poissonDefinition=DEFAULT, temperature=GRADIENT, 
        useDensity=OFF)
    compositeLayup.ReferenceOrientation(orientationType=GLOBAL, localCsys=None, 
        fieldName='', additionalRotationType=ROTATION_NONE, angle=0.0, 
        axis=AXIS_3)
    for layerNumber in range(1, len(partGrid)+1):
        layerObj = partGrid[layerNumber][objID]
        objType = layerObj.objectType
        objAngle = layerObj.angle[-1]
        if objType == 'Undulation':
            interfaceAngle = abs(layerObj.angle[0] - layerObj.angle[1])
            secAngle = (0, interfaceAngle)
            material = 'Undulation {}'.format(secAngle)
        else:
            material = objType
        nameOfPly = 'Ply-{}'.format(layerNumber)
        compositeLayup.CompositePly(suppressed=False, plyName=nameOfPly, 
            region=region, material=material, thicknessType=SPECIFY_THICKNESS, 
            thickness=0.2, orientationType=SPECIFY_ORIENT,
            orientationValue=objAngle, additionalRotationType=ROTATION_NONE,
            additionalRotationField='', axis=AXIS_3, angle=0.0, numIntPoints=3)

# create the part instance
specimenInstance = laminateAssembly.Instance(
        name='Specimen', part=specimenPart, dependent=ON)


# IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
# Mesh 

specimenFaces = specimenPart.faces
allFacesRegion = regionToolset.Region(faces=specimenFaces)
# elemType1 = mesh.ElemType(elemCode=CPE4, elemLibrary=STANDARD)
# elemType2 = mesh.ElemType(elemCode=CPE3, elemLibrary=STANDARD)
# specimenPart.setElementType(regions=allFacesRegion,
#                             elemTypes=(elemType1, elemType2))
specimenPart.setMeshControls(regions=specimenFaces, elemShape=QUAD,
                             technique=STRUCTURED)
specimenPart.seedPart(size=1.2, deviationFactor=0.1, minSizeFactor=0.1)
specimenPart.generateMesh()

# IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
# Apply boundary conditions

strainVector = [1.0, 0.0, 0.0]  # e11, e22, e12
periodicBC(modelName=modelName, xmin=xMin, ymin=yMin, xmax=xMax, ymax=yMax,
           strain=strainVector)


# IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
# Job 

mdb.Job(atTime=None, contactPrint=OFF, description='', echoPrint=OFF, 
    explicitPrecision=SINGLE, getMemoryFromAnalysis=True, historyPrint=OFF, 
    memory=90, memoryUnits=PERCENTAGE, model=modelName, modelPrint=OFF, 
    multiprocessingMode=DEFAULT, name=modelName, nodalOutputPrecision=SINGLE, 
    numCpus=1, queue=None, scratch='', type=ANALYSIS, userSubroutine='', 
    waitHours=0, waitMinutes=0)

# # Change to assembly view and hide all datums
# session.viewports['Viewport: 1'].assemblyDisplay.setValues(
#         optimizationTasks=OFF, geometricRestrictions=OFF, stopConditions=OFF)
# leaf = dgm.LeafFromDatums(laminateAssembly.datums.values())
# session.viewports['Viewport: 1'].partDisplay.displayGroup.remove(leaf=leaf)
# session.viewports['Viewport: 1'].view.setValues(session.views['Iso'])