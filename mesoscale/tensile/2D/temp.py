from sys import path
path.append('C:\\Python27\\Lib\\site-packages')
path.append('C:\\GitHub\interlaced_model_creation\\mesoscale')
from shapely.geometry import Point, Polygon, LineString
from shapely import affinity
from abaqus import *
from abaqusConstants import *
import displayGroupMdbToolset as dgm
import regionToolset
import sigc as sigc
import tapePlacement as tp
import interlacedMaterials as mats
import math


class TensileModel():
    def __init__(self, modelName):
        self.modelName = modelName
        mdb.Model(name=self.modelName, modelType=STANDARD_EXPLICIT)
        self.model = mdb.models[self.modelName]
        self.assembly = self.model.rootAssembly
        # set the work directory
        wd = 'C:\\Workspace\\{}'.format(self.modelName)
        if not os.path.exists(wd):
            os.makedirs(wd)
        os.chdir(wd)

    def setTestParameters(self, time, outputIntervals, testSpeed):
        '''Set drop weight tower test parameters'''
        self.time = time
        self.outputIntervals = outputIntervals
        self.testSpeed = testSpeed

    def setSpecimenParameters(self, t_angles, t_widths, t_spacing, t_thickness,
                              u_ratio, nPlies):
        '''
        Set interlaced specimen parameters. Note t_xxx indicates a tape
        property, u_xxx indicates an undulation property, and l_xxx indicates
        a laminate property
        '''
        self.t_angles = t_angles
        self.t_widths = [t_widths] * len(self.t_angles)
        self.t_spacing = t_spacing
        self.t_thickness = t_thickness
        self.u_ratio = u_ratio
        self.nPlies = nPlies
        self.sequences = (nPlies / int(len(self.t_angles))) / 2
        self.u_width = (self.t_thickness / self.u_ratio) / 2.0
        self.l_thickness = nPlies * self.t_thickness

    def createMaterials(self):
        '''Create material data (uses imported matData module)'''
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

        mats.tapeDamage(self.model, E11, E22, E33, nu12, nu13, nu23, G12, G13,
                        G23, density, Xt, Xc, Yt, Yc, Sl, alpha0, G1Plus,
                        G1Minus, G2Plus, G2Minus, G6)
        mats.undulationDamage(self.model, E11, E22, E33, nu12, nu13, nu23, G12,
                              G13, G23, density, Xt, Xc, Yt, Yc, Sl, alpha0,
                              G1Plus, G1Minus, G2Plus, G2Minus, G6,
                              self.t_angles, self.t_thickness, self.u_width)
        mats.undulationDamageResin(self.model, E11, E22, E33, nu12, nu13, nu23,
                                   G12, G13, G23, density, Xt, Xc, Yt, Yc, Sl,
                                   alpha0, G1Plus, G1Minus, G2Plus, G2Minus,
                                   G6, self.t_angles, self.t_thickness,
                                   self.u_width)
        mats.resinElastic(self.model)

    def createTapePaths(self, specimenDimensions):
        # create grids using shapely
        tapePaths, partGrid = tp.laminateCreation(
            tapeAngles=self.t_angles, tapeWidths=self.t_widths,
            tapeSpacing=self.t_spacing, specimenSize=specimenDimensions,
            undulationWidth=self.u_width, xMax=xMax, yMax=yMax)
        return tapePaths, partGrid

    def createPart(self):
        specimenSketch = self.model.ConstrainedSketch(name='Specimen Sketch',
                                                         sheetSize=200.0)
        specimenSketch.rectangle(point1=(xMin, yMin), point2=(xMax, yMax))
        self.specimenPart = self.model.Part(name='Specimen',
                                          dimensionality=THREE_D,
                                          type=DEFORMABLE_BODY)
        self.specimenPart.BaseShell(sketch=specimenSketch)
        self.assembly.Instance(name='Specimen Instance', part=self.specimenPart,
            dependent=ON)

    def partitionSpecimen(self):
        '''
        Function to partition the specimen into a grid of polygons depending on the
        tape angles of the laminate.

        '''
        # define mirror point (used to mirror the Polygon boundaries)
        mirrorPoint = Point([(0.0, 0.0), (0.0, 0.0)])

        tapes = zip(self.t_angles, self.t_widths)
        partitionLines = []
        uw = self.u_width
        for a, w in tapes:
            if a == 90:
                offset = w
                maxOffset = 50.0-w/2.0
                numberOffsets = int(maxOffset/offset)
                offsetList = [offset*k for k
                            in range(0, numberOffsets + 20)]
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
                offsetList = [offset*k for k in range(0, numberOffsets + 20)]
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
                dPointIDs.append(self.specimenPart.DatumPointByCoordinate(
                    coords=coordinates).id)
            dPoint1 = self.specimenPart.datums[dPointIDs[0]]
            dPoint2 = self.specimenPart.datums[dPointIDs[1]]
            dPoint3 = self.specimenPart.datums[dPointIDs[2]]
            dpbtp = self.specimenPart.DatumPlaneByThreePoints(
                point1=dPoint1, point2=dPoint2, point3=dPoint3).id
            specimenFaces = self.specimenPart.faces
            try:
                self.specimenPart.PartitionFaceByDatumPlane(
                    datumPlane=self.specimenPart.datums[dpbtp], faces=specimenFaces)
            except:
                continue
        

    def assignProperties(self, partGrid):
        for objID, obj in partGrid[1].iteritems():
            objCentroid = obj.centroid.coords
            selectedFace = self.specimenPart.faces.findAt(((objCentroid[0][0],
                                                objCentroid[0][1], 0.0), ))
            region = regionToolset.Region(faces=selectedFace)
            compositeLayup = self.specimenPart.CompositeLayup(
                name='CompositeLayup-{}'.format(objID), description='',
                elementType=SHELL, offsetType=MIDDLE_SURFACE, symmetric=False, 
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
                    print layerObj.angle
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

    def createExplicitStep(self):
        # mass scaling applied if increment < 1e-6
        self.model.ExplicitDynamicsStep(
            name='Loading Step', previous='Initial', timePeriod=self.time,
            massScaling=((SEMI_AUTOMATIC, MODEL, THROUGHOUT_STEP, 0.0, 1e-06,
                          BELOW_MIN, 1, 0, 0.0, 0.0, 0, None), ))
        self.model.fieldOutputRequests['F-Output-1'].setValues(
            variables=('S', 'LE', 'U', 'V', 'A', 'RF', 'CSTRESS', 'CSDMG',
                       'SDV', 'STATUS'),
            numIntervals=outputIntervals)

    def getFacesByBBox(self, boundingBox):
        x1, y1, z1, x2, y2, z2 = boundingBox
        faces = [inst.faces.getByBoundingBox(x1, y1, z1, x2, y2, z2)
                 for inst in self.assembly.instances.values()
                 if inst.faces.getByBoundingBox(x1, y1, z1, x2, y2, z2)]
        return faces

    def applyContraints(self, dimensions):
        tol = 0.001
        zMax = self.l_thickness
        zMin = 0.0
        xMin, yMin, xMax, yMax = dimensions
        # identify faces at top and bottom for coupling
        tFacesBBox = (xMin - tol, yMax - tol, zMin - tol, xMax + tol,
                      yMax + tol, zMax + tol)
        tFaces = self.getFacesByBBox(tFacesBBox)
        tFacesSurface = self.assembly.Surface(side1Faces=tFaces,
                                              name='Top Faces')
        bFacesBBox = (xMin - tol, yMin - tol, zMin - tol, xMax + tol,
                      yMin + tol, zMax + tol)
        bFaces = self.getFacesByBBox(bFacesBBox)
        bFacesSurface = self.assembly.Surface(side1Faces=bFaces,
                                              name='Bottom Faces')
        symFacesBBox = (xMin - tol, yMin - tol, -tol, xMax + tol,
                        yMax + tol, tol)
        symFaces = self.getFacesByBBox(symFacesBBox)
        symFacesSet = self.assembly.Set(faces=symFaces, name='Symmetry Faces')

        # create reference points for coupling
        rfPoint1Id = self.assembly.ReferencePoint(
            point=(0.0, yMax + 5.0, 0.0)).id
        rfPoint1 = self.assembly.referencePoints[rfPoint1Id]
        rfPoint1Region = self.assembly.Set(referencePoints=(rfPoint1,),
                                           name='Coupling Reference Point 1')
        rfPoint2Id = self.assembly.ReferencePoint(
            point=(0.0, yMin - 5.0, 0.0)).id
        rfPoint2 = self.assembly.referencePoints[rfPoint2Id]
        rfPoint2Region = self.assembly.Set(referencePoints=(rfPoint2,),
                                           name='Coupling Reference Point 2')

        self.model.Coupling(
            name='Top Coupling', controlPoint=rfPoint1Region,
            surface=tFacesSurface, influenceRadius=WHOLE_SURFACE,
            couplingType=KINEMATIC, localCsys=None, u1=ON, u2=ON, u3=ON,
            ur1=ON, ur2=ON, ur3=ON)

        self.model.Coupling(
            name='Bottom Coupling', controlPoint=rfPoint2Region,
            surface=bFacesSurface, influenceRadius=WHOLE_SURFACE,
            couplingType=KINEMATIC, localCsys=None, u1=ON, u2=ON,
            u3=ON, ur1=ON, ur2=ON, ur3=ON)

        self.model.SmoothStepAmplitude(name='Smoothing Amplitude', timeSpan=STEP,
                                       data=((0.0, 0.0), (1e-05, 1.0)))
        self.model.DisplacementBC(
            name='Bottom Surface BC', createStepName='Loading Step',
            region=rfPoint2Region, u1=0.0, u2=0.0, u3=0.0, ur1=0.0, ur2=0.0,
            ur3=0.0, amplitude=UNSET, fixed=OFF, distributionType=UNIFORM,
            fieldName='', localCsys=None)
        self.model.VelocityBC(
            name='Top Surface BC', createStepName='Loading Step',
            region=rfPoint1Region, v1=0.0, v2=self.testSpeed, v3=0.0,
            vr1=0.0, vr2=0.0, vr3=0.0, amplitude='Smoothing Amplitude',
            localCsys=None, distributionType=UNIFORM, fieldName='')
        self.model.ZsymmBC(name='Symmetry', createStepName='Loading Step',
                        region=symFacesSet, localCsys=None)

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

    def saveModel(self):
        mdb.saveAs(pathName=self.modelName)

    def inputPopup(self):
        userInput = getInputs(fields=(('Laminate Angles (deg):', '(0, 90)'),
                                      ('Tape Width (mm):', '25'),
                                      ('Tape Spacing (int):', '1'),
                                      ('Cured Ply Thickness (mm):', '0.18'),
                                      ('Undulation Ratio (float):', '0.09')),
                              label='Please provide the following information',
                              dialogTitle='Model Parameters')
        self.t_angles = ast.literal_eval(userInput[0])  # angles of tapes
        tapeWidth = float(userInput[1])
        self.t_widths = (tapeWidth, ) * len(self.t_angles)  # tape widths
        # number of gaps between tapes in interlacing pattern
        self.t_spacing = int(userInput[2])
        # cured ply thickness e.g. 0.18125
        self.t_thickness = float(userInput[3])
        # ratio of undulation amplitude to length e.g. 0.090625
        self.u_ratio = float(userInput[4])

if __name__ == '__main__':

    # Simulation parameters
    time = 5.0  # duration to simulate [ms]
    outputIntervals = 50  # requested field output intervals
    testSpeed = 0.5  # crosshead velocity

    # Material parameters
    tapeAngles = (0, 90)  # define angles of tapes
    tapeWidths = 15.0
    tapeSpacing = 1  # number of gaps between tapes in interlacing pattern
    tapeThickness = 0.18  # cured ply thickness e.g. 0.18125
    undulationRatio = 0.09  # ratio of undulation amplitude to length
    nLayers = 4  # symmetric 8 ply laminate

    # Mesh sizes
    fineMesh = 0.125
    mediumMesh = 1.0
    coarseMesh = 3.0

    # Define specimen dimensions
    xMin = -25.0
    xMax = -xMin
    yMin = -75.0
    yMax = -yMin
    fineRegion = Polygon([(xMin, yMin), (xMax, yMin), (xMax, yMax),
                          (xMin, yMax)])
    dimensions = [xMin, yMin, xMax, yMax]

    # RVE dimensions
    # xMin = yMin = -(tapeWidths / 2.0)
    # xMax = yMax = xMin + (tapeSpacing + 1) * (tapeWidths)
    # yMin = -75.0
    # yMax = 75.0
    # fineRegion = Polygon([(xMin, yMin), (xMax, yMin), (xMax, yMax),
    #                      (xMin, yMax)])
    # dimensions = [xMin, yMin, xMax, yMax]
    

    # Create model
    mdl = TensileModel('Test')
    mdl.setTestParameters(time, outputIntervals, testSpeed)
    mdl.setSpecimenParameters(tapeAngles, tapeWidths, tapeSpacing,
                              tapeThickness, undulationRatio, nLayers)
    mdl.createMaterials()
    paths_f, grid_f = mdl.createTapePaths(fineRegion)
    mdl.createPart()
    mdl.partitionSpecimen()
    mdl.assignProperties(grid_f)
    mdl.createExplicitStep()
    mdl.applyContraints(dimensions)
    mdl.createJob()
    mdl.saveModel()
