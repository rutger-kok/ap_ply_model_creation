'''
Module to determine the warpage of an interlaced laminate using Abaqus.
2D shell model with a predefined -150 degree temperature change. The
interlaced laminate does not include undulation regions.

Not yet updated to conform with class inheritance strategy.

(c) Rutger Kok 03/12/2020
'''

from abaqus import *
from abaqusConstants import *
from caeModules import *
from visualization import *
import mesh
import regionToolset
import os
from sys import path
path.append('C:\\Python27\\Lib\\site-packages')
path.append('C:\\Workspace\\3D_Mechprops')
from shapely.geometry import Polygon, Point, LineString
from shapely.affinity import scale, rotate
from shapely.ops import linemerge, unary_union, polygonize
from math import cos, radians, tan, sin, pi
import interlacedMaterials2 as mats


def createGrids(tapeAngles, tapeWidths, sample=None):
    '''
    Function to partition the specimen into a grid of polygons depending on the
    tape angles of the laminate.
    The specimen region (150mm x 100mm) is separated into Shapely Polygons.
    '''
    # define specimen boundaries
    if sample is None:
        sample = Polygon(
            [(-50.0, -75.0), (50.0, -75.0), (50.0, 75.0), (-50.0, 75.0)])

    xMax = sample.exterior.xy[0][2]
    yMax = sample.exterior.xy[1][2]

    # define mirror point (used to mirror the Polygon boundaries)
    mirrorPoint = Point([(0.0, 0.0), (0.0, 0.0)])

    tapes = zip(tapeAngles, tapeWidths)
    partitionLines = []
    for tape in tapes:
        a = tape[0]  # angle
        w = tape[1]  # width
        if a == 90:
            offset = w
            maxOffset = xMax - w / 2.0
            numberOffsets = int(maxOffset / offset)
            offsetList = [offset * k for k in range(0, numberOffsets + 1)]
            for os in offsetList:
                tapeLineCoords = [((w / 2.0) + os, 900.0),
                                  ((w / 2.0) + os, -900.0)]
                tapeLine = LineString(tapeLineCoords)
                partitionLines.append(tapeLine)
            reflectedLines = [scale(line, xfact=-1, origin=mirrorPoint)
                              for line in partitionLines]
            partitionLines = partitionLines + reflectedLines

        else:
            offset = w / cos(radians(a))
            maxOffset = (yMax - (w / 2.0) * cos(radians(a))
                         + xMax * tan(radians(a)))
            numberOffsets = int(maxOffset / offset)
            offsetList = [offset * k for k in range(0, numberOffsets + 20)]
            for os in offsetList:
                tapeLineCoords = [(-900.0, (w / 2.0) + os),
                                  (900.0, (w / 2.0) + os)]
                rotPoint = Point([(0.0, os), (0.0, os)])
                tapeLine = rotate(LineString(tapeLineCoords), a, rotPoint)
                partitionLines.append(tapeLine)
            reflectedLines = [rotate(line, 180.0, origin=mirrorPoint)
                              for line in partitionLines]
            partitionLines = partitionLines + reflectedLines

    # collection of individual linestrings for splitting in a list and add
    # the polygon lines to it.
    partitionLines.append(sample.boundary)
    merged_lines = linemerge(partitionLines)
    border_lines = unary_union(merged_lines)
    decomposition = polygonize(border_lines)
    trimmed = [plygn for plygn in decomposition if plygn.within(sample)]
    polyDict = {key: value for (key, value) in enumerate(trimmed)}
    return polyDict


# initialize additional polygon attributes
setattr(Polygon, 'angle', None)  # create attribute in Polygon class
setattr(Polygon, 'A', None)  # create attribute in Polygon class
setattr(Polygon, 'B', None)  # create attribute in Polygon class
setattr(Polygon, 'D', None)  # create attribute in Polygon class


def machinePass(grid, angle=0, coords=None):
    '''Function defining a single pass of a an AFP machine
    inputs: grid of the specimen, angle of the pass, coords
    NOTE: a pass should be specified using its coordinates when horizontal
    and then rotated. Do not rotate outside of the function'''

    # default tape dimensions (used unless other coords are supplied)
    if coords is None:
        w = 25.0
        coords = ([(-300.0, (w / 2)), (300.0, (w / 2)),
                   (300.0, (-w / 2)), (-300.0, (-w / 2))])
    else:
        coords = coords

    # define tape boundaries
    bounds = Polygon(coords)
    rotatedBounds = rotate(bounds, angle).buffer(1 * 10**-8, join_style=2)
    # identify which Polygons from the createGrids function are within
    # the boundaries of the rotated tape
    tapePath = [(obj1ID, obj1) for (obj1ID, obj1) in grid.iteritems()
                if obj1.within(rotatedBounds)]

    for i in range(len(tapePath)):
        obj2ID, obj2 = tapePath[i]
        if obj2.angle is None:
            obj2.angle = [angle]
        elif isinstance(obj2.angle, list):
            obj2.angle.append(angle)
        else:
            print 'Error: obj.angle type is = '.format(type(obj2.angle))


def laminateCreation(grid, tapeAngles, tapeWidths, tapeSpacing):
    tapes = zip(tapeAngles, tapeWidths)
    s = tapeSpacing  # spacing
    for numShift in range(s + 1):
        for a, w in tapes:
            if a == 90:
                offset = w
                spacedOffset = (1 + s) * offset
                maxOffset = 50.0 + w / 2.0
                numberOffsets = int(maxOffset / spacedOffset)
                offsetList = [spacedOffset * k for k
                              in range(-numberOffsets - numShift - 10,
                                       numberOffsets - numShift + 10)]
                shift = numShift * offset
                for os in offsetList:
                    tapeCoords = ([(-400.0 + os + shift, (w / 2.0)),
                                   (400.0 + os + shift, (w / 2.0)),
                                   (400.0 + os + shift, (-w / 2.0)),
                                   (-400.0 + os + shift, (-w / 2.0))])
                    machinePass(grid, coords=tapeCoords, angle=90)
            else:
                offset = w / cos(radians(a))
                spacedOffset = ((1 + s) * w) / cos(radians(a))
                maxOffset = ((150.0 + w / 2.0) * cos(radians(a))
                             + 150.0 * tan(radians(a)))
                numberOffsets = int(maxOffset / spacedOffset)
                offsetList = [spacedOffset * k for k
                              in range(-numberOffsets - numShift - 10,
                                       numberOffsets - numShift + 10)]
                shift = numShift*offset
                for os in offsetList:
                    tapeCoords = ([(-400.0, ((w / 2.0)) + os + shift),
                                   (400.0, ((w / 2.0)) + os + shift),
                                   (400.0, ((-w / 2.0)) + os + shift),
                                   (-400.0, ((-w / 2.0)) + os + shift)])
                    machinePass(grid, coords=tapeCoords, angle=a)


class TensileModel():
    def __init__(self):
        self.modelName = 'QI_Interlace'
        mdb.Model(name=self.modelName, modelType=STANDARD_EXPLICIT)
        self.model = mdb.models[self.modelName]
        self.assembly = self.model.rootAssembly
        # set the work directory
        wd = 'C:\\Workspace\\Warpage\\{}'.format(self.modelName)
        if not os.path.exists(wd):
            os.makedirs(wd)
        os.chdir(wd)


    def setSpecimenParameters(self, t_angles, t_widths, t_spacing, t_thickness,
                              nPlies):
        '''
        Set interlaced specimen parameters. Note t_xxx indicates a tape
        property, u_xxx indicates an undulation property, and l_xxx indicates
        a laminate property
        '''
        self.t_angles = t_angles
        self.t_widths = [t_widths] * len(self.t_angles)
        self.t_spacing = t_spacing
        self.t_thickness = t_thickness
        self.nPlies = nPlies
        self.sequences = (nPlies / int(len(self.t_angles)))
        self.l_thickness = nPlies * self.t_thickness

    def createMaterials(self):
        '''Create material data (uses imported module)'''
        E11 = 146.8
        E22 = E33 = 11.6
        nu12 = nu13 = 0.289
        nu23 = 0.298
        G12 = G13 = 6.47
        G23 = 4.38
        Xt = 2.354
        Xc = 1.102
        Yt = 0.0343
        Yc = 0.184
        Sl = 0.0827
        alpha0 = 53.0
        G1Plus = 0.1
        G1Minus = 0.1
        G2Plus = 0.00075
        G2Minus = 0.0025
        G6 = 0.0035
        density = 1.59e-06

        mats.tapeElastic(self.model, E11, E22, E33, nu12, nu13, nu23, G12, G13,
                         G23, density)
        self.model.materials['Tape-Elastic'].Expansion(type=ORTHOTROPIC, 
            table=((-9e-7, 2.2e-5, 2.2e-5), ))


    def createTapePaths(self, specimenSize):
        # create grids using shapely
        grid = createGrids(tapeAngles=self.t_angles, tapeWidths=self.t_widths, sample=specimenSize)
        tapePaths = laminateCreation(grid, tapeAngles=self.t_angles,
                                    tapeWidths=self.t_widths, tapeSpacing=self.t_spacing)
        return tapePaths, grid

    def partitionSpecimen(self):
        '''
        Function to partition the specimen into a grid of polygons depending on the
        tape angles of the laminate.

        '''
        # define mirror point (used to mirror the Polygon boundaries)
        mirrorPoint = Point([(0.0, 0.0), (0.0, 0.0)])

        partitionLines = []
        for a, w in zip(self.t_angles, self.t_widths):
            if a == 90:
                offset = w
                maxOffset = xMax - w / 2.0
                numberOffsets = int(maxOffset / offset)
                offsetList = [offset * k for k in range(0, numberOffsets + 1)]
                for os in offsetList:
                    tapeLineCoords = [((w / 2.0) + os, 300.0),
                                    ((w / 2.0) + os, -300.0)]
                    tapeLine = LineString(tapeLineCoords)
                    partitionLines.append(tapeLine)
                reflectedLines = [scale(line, xfact=-1, origin=mirrorPoint)
                                for line in partitionLines]
                partitionLines = partitionLines + reflectedLines

            else:
                offset = w / cos(radians(a))
                maxOffset = (yMax - (w / 2.0) * cos(radians(a))
                            + xMax * tan(radians(a)))
                numberOffsets = int(maxOffset / offset)
                offsetList = [offset * k for k in range(0, numberOffsets + 20)]
                for os in offsetList:
                    tapeLineCoords = [(-300.0, (w / 2.0) + os),
                                    (300.0, (w / 2.0) + os)]
                    rotPoint = Point([(0.0, os), (0.0, os)])
                    tapeLine = rotate(LineString(tapeLineCoords), a, rotPoint)
                    partitionLines.append(tapeLine)
                reflectedLines = [rotate(line, 180.0, origin=mirrorPoint)
                                for line in partitionLines]
                partitionLines = partitionLines + reflectedLines

        specimenPart = self.specimenPart
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

    def createPart(self):
        specimenSketch = self.model.ConstrainedSketch(name='Specimen Sketch',
                                                         sheetSize=200.0)
        specimenSketch.rectangle(point1=(xMin, yMin), point2=(xMax, yMax))
        self.specimenPart = self.model.Part(name='Specimen',
                                          dimensionality=THREE_D,
                                          type=DEFORMABLE_BODY)
        self.specimenPart.BaseShell(sketch=specimenSketch)
        self.specimenFaces = self.specimenPart.faces

        faceRegion = regionToolset.Region(faces=self.specimenFaces[:])
        self.assembly.Instance(name='Specimen Instance', part=self.specimenPart,
            dependent=ON)

    def assignProperties(self, partGrid):
        angleSets = set()
        for objID, obj in partGrid.iteritems():
            angleSets.add(tuple(obj.angle))

        for k, angle in enumerate(angleSets):
            centroids = [obj.centroid.coords for obj in partGrid.values() if tuple(obj.angle) == angle]
            selectedFaces = self.specimenFaces.findAt(((centroids[0][0][0], centroids[0][0][1], 0.0), ))
            for n in range(1, len(centroids)):
                selectedFaces += self.specimenFaces.findAt(((centroids[n][0][0], centroids[n][0][1], 0.0), ))
            region = regionToolset.Region(faces=selectedFaces)
            compositeLayup = self.specimenPart.CompositeLayup(
                name='CompositeLayup-{}'.format(k), description='',
                elementType=SHELL, offsetType=MIDDLE_SURFACE, symmetric=False,
                thicknessAssignment=FROM_SECTION)
            compositeLayup.Section(preIntegrate=OFF, integrationRule=SIMPSON,
                thicknessType=UNIFORM, poissonDefinition=DEFAULT, temperature=GRADIENT,
                useDensity=OFF)
            compositeLayup.ReferenceOrientation(orientationType=GLOBAL, localCsys=None,
                fieldName='', additionalRotationType=ROTATION_NONE, angle=0.0,
                axis=AXIS_3)
            for i in range(self.sequences):
                for j, subAngle in enumerate(angle):
                    nameOfPly = 'Ply-{}-{}'.format(i,j)
                    compositeLayup.CompositePly(suppressed=False, plyName=nameOfPly,
                        region=region, material='Tape-Elastic', thicknessType=SPECIFY_THICKNESS,
                        thickness=self.t_thickness, orientationType=SPECIFY_ORIENT,
                        orientationValue=subAngle, additionalRotationType=ROTATION_NONE,
                        additionalRotationField='', axis=AXIS_3, angle=0.0, numIntPoints=3)

        # partition part to facilitate central BC
        dplaneID1 = self.specimenPart.DatumPlaneByPrincipalPlane(
            principalPlane=YZPLANE, offset=(0.0)).id
        self.specimenPart.PartitionFaceByDatumPlane(
            datumPlane=self.specimenPart.datums[dplaneID1],
            faces=self.specimenPart.faces)
        dplaneID2 = self.specimenPart.DatumPlaneByPrincipalPlane(
            principalPlane=XZPLANE, offset=(0.0)).id
        self.specimenPart.PartitionFaceByDatumPlane(
            datumPlane=self.specimenPart.datums[dplaneID2],
            faces=self.specimenPart.faces)

        # mesh the part
        faceRegion = regionToolset.Region(faces=self.specimenFaces[:])
        elemType1 = mesh.ElemType(
            elemCode=S4R, elemLibrary=STANDARD, secondOrderAccuracy=OFF,
            hourglassControl=ENHANCED)
        elemType2 = mesh.ElemType(elemCode=S3, elemLibrary=STANDARD)
        self.specimenPart.setElementType(regions=faceRegion, elemTypes=(elemType1, elemType2))
        self.specimenPart.seedPart(size=2.5)
        self.specimenPart.generateMesh()


    def createImplicitStep(self):
        self.model.StaticStep(name='Loading Step', previous='Initial')


    def applyContraints(self, dimensions):

        # all faces
        specimenInstance = self.assembly.instances['Specimen Instance']
        allFaces = specimenInstance.faces
        allFacesSet = self.assembly.Set(faces=allFaces, name='All Faces')
        self.model.Temperature(
            name='HighTemp', createStepName='Initial', region=allFacesSet,
            distributionType=UNIFORM, magnitudes=(-150.0, ),
            crossSectionDistribution=CONSTANT_THROUGH_THICKNESS)

        self.model.Temperature(
            name='LowTemp', createStepName='Loading Step', region=allFacesSet,
            distributionType=UNIFORM, magnitudes=(0.0, ),
            crossSectionDistribution=CONSTANT_THROUGH_THICKNESS)

        fixPoint = specimenInstance.vertices.findAt(((0.0, 0.0, 0.0), ))
        fixPointRegion = self.assembly.Set(vertices=fixPoint,
                                           name='Fixed Point')

        self.model.DisplacementBC(
            name='Fix Point', createStepName='Loading Step',
            region=fixPointRegion, u1=0.0, u2=0.0, u3=0.0, ur1=UNSET, ur2=UNSET,
            ur3=UNSET, amplitude=UNSET, fixed=OFF, distributionType=UNIFORM,
            fieldName='', localCsys=None)


    def createJob(self, cpus=4):
        mdb.Job(
            name=self.modelName, model=self.modelName, description='',
            atTime=None, waitMinutes=0, waitHours=0, queue=None, type=ANALYSIS,
            memoryUnits=PERCENTAGE, explicitPrecision=DOUBLE_PLUS_PACK,
            nodalOutputPrecision=FULL, echoPrint=OFF, modelPrint=OFF,
            contactPrint=OFF, historyPrint=OFF, userSubroutine='', scratch='',
            resultsFormat=ODB, parallelizationMethodExplicit=DOMAIN,
            numDomains=cpus, activateLoadBalancing=False,
            multiprocessingMode=DEFAULT, numCpus=cpus, memory=90)
        # mdb.jobs[self.modelName].submit(consistencyChecking=OFF)
        # mdb.jobs[self.modelName].waitForCompletion()

    def saveModel(self):
        mdb.saveAs(pathName=self.modelName)


if __name__ == '__main__':

    # Material parameters
    tapeAngles = (0, 45, 90, -45)  # define angles of tapes
    tapeWidths = 10.0
    tapeSpacing = 3  # number of gaps between tapes in interlacing pattern
    tapeThickness = 0.2  # cured ply thickness e.g. 0.18125
    nLayers = 28

    # Mesh sizes
    fineMesh = 0.125
    mediumMesh = 1.0
    coarseMesh = 3.0

    # Define specimen dimensions
    # xMin = -25.0
    # xMax = -xMin
    # yMin = -75.0
    # yMax = -yMin
    # fineRegion = Polygon([(xMin, yMin), (xMax, yMin), (xMax, yMax),
    #                       (xMin, yMax)])
    # dimensions = [xMin, yMin, xMax, yMax]

    # RVE dimensions
    xMin = -150.0
    xMax = 150.0
    yMin = -150.0
    yMax = 150.0
    fineRegion = Polygon([(xMin, yMin), (xMax, yMin), (xMax, yMax),
                         (xMin, yMax)])
    dimensions = [xMin, yMin, xMax, yMax]
    

    # Create model
    mdl = TensileModel()
    mdl.setSpecimenParameters(tapeAngles, tapeWidths, tapeSpacing,
                              tapeThickness, nLayers)
    mdl.createMaterials()
    paths_f, grid_f = mdl.createTapePaths(fineRegion)
    mdl.createPart()
    mdl.partitionSpecimen()
    mdl.assignProperties(grid_f)
    mdl.createImplicitStep()
    mdl.applyContraints(dimensions)
    mdl.createJob()
    mdl.saveModel()
