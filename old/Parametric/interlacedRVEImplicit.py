from sys import path
githubPath = r"\\arran.sms.ed.ac.uk\home\s1342398\GitHub"
path.append('C:\\Python27\\Lib\\site-packages')
path.append(githubPath + '\\interlaced_model_creation\\editing')
path.append(githubPath + '\\interlaced_model_creation\\editing\\3D\\Mechanical Properties\\Parametric')
from shapely.geometry import Point, Polygon, LineString
from shapely import affinity
from abaqus import *
from abaqusConstants import *
from visualization import *
import regionToolset
import sigc as sigc
import tapePlacement as tp
from math import tan, cos, radians
from periodicBC_3D import periodicBC
import os
import interlacedMaterialProps as matprops
from caeModules import *

def createPartGrids(tapeAngles, tapeWidths, tapeSpacing, tapeThickness,
                    undulationWidth=1.0):
    # RVE dimensions
    xMin = yMin = -(tapeWidths[0] / 2.0)
    xMax = yMax = xMin + (tapeSpacing + 1) * (tapeWidths[0])
    numLayers = len(tapeAngles) * 2.0  # symmetric
    laminateThickness = numLayers * tapeThickness
    zMax = laminateThickness / 2.0
    zMin = -zMax
    dimensions = [xMin, yMin, zMin, xMax, yMax, zMax]

    RVEPolygon = Polygon([(xMin, yMin), (xMax, yMin), (xMax, yMax),
                          (xMin, yMax)])

    # create grid using shapely
    partGrid = sigc.createGrids(tapeAngles=tapeAngles, tapeWidths=tapeWidths,
                                undulationWidth=undulationWidth,
                                sample=RVEPolygon)
    # identify parts in grid
    tapePaths = tp.laminateCreation(grid=partGrid, tapeAngles=tapeAngles,
                                    tapeWidths=tapeWidths,
                                    tapeSpacing=tapeSpacing,
                                    xMax=xMax, yMax=yMax,
                                    undulationWidth=undulationWidth)

    return partGrid, tapePaths, dimensions


def createStaticStep(modelName):
    activeModel = mdb.models[modelName]
    activeModel.StaticStep(name='Loading Step', previous='Initial',
                           initialInc=0.1)


def createSpecimenPart(modelName, partName, dimensions):
    activeModel = mdb.models[modelName]

    p1 = (dimensions[0], dimensions[1])
    p2 = (dimensions[3], dimensions[4])
    d = dimensions[5]

    # Create sketch
    partSketch = activeModel.ConstrainedSketch(name='Part Sketch',
                                               sheetSize=200.0)
    partSketch.rectangle(point1=p1, point2=p2)
    part = activeModel.Part(name=partName, dimensionality=THREE_D,
                            type=DEFORMABLE_BODY)
    part.BaseSolidExtrude(sketch=partSketch, depth=d)


def partitionPart(modelName, partName, dimensions, tapeAngles, tapeWidths,
                  tapeThickness, undulationWidth=1.0):
    '''
    Function to partition the specimen into plies, and subsequently into a grid
    of polygons depending on the tape angles of the laminate.
    '''
    activeModel = mdb.models[modelName]
    activePart = activeModel.parts[partName]
    xMax = dimensions[3]
    yMax = dimensions[4]

    # Through thickness partition into plies
    plyDatumPlaneIds = []
    for k in range(1, len(tapeAngles)):
        plyDatumPlaneIds.append(activePart.DatumPlaneByPrincipalPlane(
            principalPlane=XYPLANE, offset=tapeThickness * k).id)
        partCells = activePart.cells
        activePart.PartitionCellByDatumPlane(
            datumPlane=activePart.datums[plyDatumPlaneIds[k - 1]],
            cells=partCells)

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
            offset = w / cos(radians(a))
            maxOffset = (yMax - (w / 2.0) * cos(radians(a)) +
                         50.0 * tan(radians(a)))
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
            dPointIDs.append(activePart.DatumPointByCoordinate(
                coords=coordinates).id)
        dPoint1 = activePart.datums[dPointIDs[0]]
        dPoint2 = activePart.datums[dPointIDs[1]]
        dPoint3 = activePart.datums[dPointIDs[2]]
        dpbtp = activePart.DatumPlaneByThreePoints(point1=dPoint1,
                                                   point2=dPoint2,
                                                   point3=dPoint3).id
        specimenCells = activePart.cells
        try:
            activePart.PartitionCellByDatumPlane(
                datumPlane=activePart.datums[dpbtp], cells=specimenCells)
        except:
            continue


def mirrorPart(modelName, partName):
    activeModel = mdb.models[modelName]
    activePart = activeModel.parts[partName]
    mirrorFace = activePart.faces.findAt(coordinates=(0.0, 0.0, 0.0))
    activePart.Mirror(mirrorPlane=mirrorFace, keepOriginal=ON,
                      keepInternalBoundaries=ON)


def createMesh(modelName, partName, meshSize):
    activeModel = mdb.models[modelName]
    activePart = activeModel.parts[partName]
    activePart.seedPart(size=meshSize)
    activePart.generateMesh()


def assignProps(modelName, partName, partGrid, tapeThickness):
    activeModel = mdb.models[modelName]
    activePart = activeModel.parts[partName]
    partCells = activePart.cells

    for layer, polyDict in partGrid.iteritems():
        for objID, obj in polyDict.iteritems():
            objCentroid = obj.centroid.coords
            zCoord = (2 * layer - 1) / 2.0 * tapeThickness
            selectedCells = partCells.findAt(((objCentroid[0][0],
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
            activePart.SectionAssignment(region=objRegion, sectionName=secName,
                                         offset=0.0, offsetType=MIDDLE_SURFACE,
                                         offsetField='',
                                         thicknessAssignment=FROM_SECTION)
            activePart.MaterialOrientation(region=objRegion,
                                           orientationType=SYSTEM, axis=AXIS_3,
                                           localCsys=None, fieldName='',
                                           additionalRotationType=ROTATION_ANGLE,
                                           additionalRotationField='',
                                           angle=objAngle,
                                           stackDirection=STACK_3)


def createInstance(modelName, partName, instanceName=None):
    activeModel = mdb.models[modelName]
    activeAssembly = activeModel.rootAssembly
    if instanceName is None:  # instanceName defaults to partName if undefined
        instanceName = partName
    # create the part instance
    partInstance = activeAssembly.Instance(name=instanceName,
                                           part=activeModel.parts[partName],
                                           dependent=ON)


def assignPBCs(modelName, dimensions, displacements):
    activeModel = mdb.models[modelName]
    activeAssembly = activeModel.rootAssembly
    periodicBC(modelName, dimensions)

    for ind, s in enumerate(displacements):
        dis = [0.0, ] * 3
        if s == UNSET:
            continue
        elif s:
            dis[ind] = s
        bcRegion = activeAssembly.sets['MasterNode{}'.format(ind + 1)]
        activeModel.DisplacementBC(name='BC-{}'.format(ind),
                                   createStepName='Loading Step',
                                   region=bcRegion, u1=dis[0], u2=dis[1],
                                   u3=dis[2], ur1=UNSET, ur2=UNSET,
                                   ur3=UNSET, amplitude=UNSET,
                                   localCsys=None, distributionType=UNIFORM,
                                   fieldName='')


def createJob(modelName, run=False, cpus=6):

    mdb.Job(name=modelName, model=modelName, description='', type=ANALYSIS, 
            atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90,
            memoryUnits=PERCENTAGE, getMemoryFromAnalysis=True, 
            explicitPrecision=SINGLE, nodalOutputPrecision=FULL,
            echoPrint=OFF, modelPrint=OFF, contactPrint=OFF, historyPrint=OFF,
            userSubroutine='', scratch='', resultsFormat=ODB,
            multiprocessingMode=DEFAULT, numCpus=cpus, numDomains=cpus,
            numGPUs=0)
    if run is True:
        mdb.jobs[modelName].submit(consistencyChecking=OFF)
        mdb.jobs[modelName].waitForCompletion()


def postProcess(modelName, dimensions):
    '''
    Determine the stiffness of the RVE
    '''
    xMin, yMin, zMin, xMax, yMax, zMax = dimensions
    specimenHeight = yMax - yMin
    specimenWidth = xMax - xMin
    laminateThickness = zMax - zMin

    # Open the Output Database for the current Job
    odb = openOdb(path='{}.odb'.format(modelName))

    frame = odb.steps['Loading Step'].frames[-1]
    rfFieldOutput = frame.fieldOutputs['RF']
    uFieldOutput = frame.fieldOutputs['U']
    masterNode1 = odb.rootAssembly.nodeSets['MASTERNODE1']
    rfMasterNode1 = rfFieldOutput.getSubset(region=masterNode1)
    uMasterNode1 = uFieldOutput.getSubset(region=masterNode1)

    rf1 = rfMasterNode1.values[0].dataDouble[0]
    u1 = uMasterNode1.values[0].dataDouble[0]

    E11 = ((rf1 / (laminateThickness * specimenHeight)) / (u1 / specimenWidth))

    return E11


def main(tapeAngles, tapeWidths, tapeSpacing, tapeThickness, undulationRatio,
         displacements):

    # Name model and change working directory
    ratioString = str(undulationRatio).replace('.', '_')
    angleString = '-'.join([str(x) for x in tapeAngles])
    modelName = '{}--{}-{}-{}'.format(angleString, int(tapeWidths[0]),
                                      int(tapeSpacing), ratioString)
    mdb.Model(name=modelName, modelType=STANDARD_EXPLICIT)

    # set the work directory
    wd = 'C:\\Workspace\\3D_RVE\\Parametric\\{}'.format(modelName)
    if not os.path.exists(wd):
        os.makedirs(wd)
    os.chdir(wd)

    # Create model
    undulationWidth = tapeThickness / undulationRatio
    partName = 'Specimen'
    partGrid, tapePaths, dimensions = createPartGrids(tapeAngles, tapeWidths,
                                                      tapeSpacing,
                                                      tapeThickness,
                                                      undulationWidth)
    createSpecimenPart(modelName, partName, dimensions)
    partitionPart(modelName, partName, dimensions, tapeAngles, tapeWidths,
                  tapeThickness, undulationWidth)
    mirrorPart(modelName, partName)
    matprops.VTC401_Elastic(modelName, tapeAngles, tapeThickness,
                            undulationWidth)
    assignProps(modelName, partName, partGrid, tapeThickness)
    createMesh(modelName, partName, 0.5)
    createStaticStep(modelName)
    createInstance(modelName, partName)
    assignPBCs(modelName, dimensions, displacements)
    createJob(modelName, run=True)
    stiffness = postProcess(modelName, dimensions)

    return stiffness


if __name__ == '__main__':

    # Test model parameters
    tapeAngles = (0, 90)  # define angles of tapes
    tapeWidth = 25.0
    tapeWidths = (tapeWidth, ) * len(tapeAngles)  # tape widths
    # number of gaps between tapes in interlacing pattern
    tapeSpacing = 2
    # cured ply thickness e.g. 0.18125
    tapeThickness = 0.18
    # ratio of undulation amplitude to length e.g. 0.090625
    undulationRatio = 0.18
    displacements = [0.1, UNSET, UNSET]

    main(tapeAngles, tapeWidths, tapeSpacing, tapeThickness, undulationRatio,
         displacements)
