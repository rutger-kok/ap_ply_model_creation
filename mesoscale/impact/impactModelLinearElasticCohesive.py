'''
Creates an Abaqus model 
'''

from abaqus import *
from abaqusConstants import *
import mesh
import regionToolset
import os
from math import cos, tan, radians
from sys import path
githubPath = 'C:\\GitHub'
path.append('C:\\Python27\\Lib\\site-packages')
path.append(githubPath + '\\interlaced_model_creation\\mesoscale')
from shapely.geometry import Point, Polygon, LineString
import sigc as sigc
from shapely.affinity import scale, rotate
import interlacedMaterials as mats
import tapePlacement as tp


class ImpactModel():
    def __init__(self, modelName):
        self.modelName = modelName
        mdb.Model(name=self.modelName, modelType=STANDARD_EXPLICIT)
        self.model = mdb.models[self.modelName]
        self.assembly = self.model.rootAssembly
        # set the work directory
        wd = 'C:\\Workspace\\Impact\\{}'.format(self.modelName)
        if not os.path.exists(wd):
            os.makedirs(wd)
        os.chdir(wd)

    def setTestParameters(self, time, outputIntervals, energy):
        '''Set drop weight tower test parameters'''
        self.time = time
        self.outputIntervals = outputIntervals
        self.energy = energy
        self.impactorMass = 5.0  # kg
        # calculate initial impactor velocity
        self.initialVelocity = -((2.0 * self.energy) / self.impactorMass)**0.5

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
        self.u_width = self.t_thickness / self.u_ratio
        self.l_thickness = nPlies * self.t_thickness

    def createMaterials(self):
        '''Create material data (uses imported matData module)'''
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

        mats.tapeElastic(self.model, E11, E22, E33, nu12, nu13, nu23, G12, G13,
                         G23, density)
        mats.tapeDamage(self.model, E11, E22, E33, nu12, nu13, nu23, G12, G13,
                        G23, density, Xt, Xc, Yt, Yc, Sl, alpha0, G1Plus,
                        G1Minus, G2Plus, G2Minus, G6)
        mats.undulationDamage(self.model, E11, E22, E33, nu12, nu13, nu23, G12,
                              G13, G23, density, Xt, Xc, Yt, Yc, Sl, alpha0,
                              G1Plus, G1Minus, G2Plus, G2Minus, G6,
                              self.t_angles, self.t_thickness, self.u_width)
        mats.undulationElastic(self.model, E11, E22, E33, nu12, nu13, nu23,
                               G12, G13, G23, density, self.t_angles,
                               self.t_thickness, self.u_width)
        mats.undulationDamageResin(self.model, E11, E22, E33, nu12, nu13, nu23,
                                   G12, G13, G23, density, Xt, Xc, Yt, Yc, Sl,
                                   alpha0, G1Plus, G1Minus, G2Plus, G2Minus,
                                   G6, self.t_angles, self.t_thickness,
                                   self.u_width)
        mats.undulationElasticResin(self.model, E11, E22, E33, nu12, nu13,
                                    nu23, G12, G13, G23, density,
                                    self.t_angles, self.t_thickness,
                                    self.u_width)
        mats.resinElastic(self.model)

    def createImpactorPart(self):
        ''' Creates the hemispherical impactor for DWT simulations'''
        # Sketch the impactor (to create part by revolution)
        sketch = self.model.ConstrainedSketch(name='Impactor Sketch',
                                              sheetSize=200.0)
        sketch.ConstructionLine(point1=(0.0, -100.0), point2=(0.0, 100.0))
        sketch.ArcByCenterEnds(center=(0.0, 0.0), point1=(8.0, 0.0),
                               point2=(0.0, -8.0), direction=CLOCKWISE)
        sketch.Line(point1=(0.0, -8.0), point2=(0.0, 32.0))
        sketch.Line(point1=(0.0, 32.0), point2=(8.0, 32.0))
        sketch.Line(point1=(8.0, 32.0), point2=(8.0, 0.0))
        self.impactorPart = self.model.Part(
            name='Impactor', dimensionality=THREE_D,
            type=DISCRETE_RIGID_SURFACE)
        self.impactorPart.BaseSolidRevolve(sketch=sketch, angle=360.0,
                                           flipRevolveDirection=OFF)
        RFPointID = self.impactorPart.ReferencePoint(point=(0.0, 0.0, 0.0)).id
        RFPoint = self.impactorPart.referencePoints[RFPointID]
        RFPointRegion = self.impactorPart.Set(
            referencePoints=(RFPoint,), name='Impactor Reference Point')
        # Create shell
        cells = self.impactorPart.cells
        self.impactorPart.RemoveCells(cellList=cells[0:1])
        # Assign inertia/mass to impactor
        self.impactorPart.engineeringFeatures.PointMassInertia(
            name='Impactor Mass-Inertia', region=RFPointRegion,
            mass=self.impactorMass, alpha=0.0, composite=0.0)
        # Meshing
        elemType1 = mesh.ElemType(elemCode=R3D4, elemLibrary=STANDARD)
        elemType2 = mesh.ElemType(elemCode=R3D3, elemLibrary=STANDARD)
        faces = self.impactorPart.faces
        pickedRegions = (faces, )
        self.impactorPart.setElementType(regions=pickedRegions,
                                         elemTypes=(elemType1, elemType2))
        self.impactorPart.seedPart(size=1.0, deviationFactor=0.1,
                                   minSizeFactor=0.1)
        self.impactorPart.generateMesh()

    def createClampParts(self):
        '''Create the restraining clamps for DWT simulations'''
        coords = (43.75, 50.0)
        sketch = self.model.ConstrainedSketch(name='Clamp Sketch',
                                              sheetSize=200.0)
        sketch.CircleByCenterPerimeter(center=coords,
                                       point1=(coords[0] + 4.0, coords[1]))
        self.clampPart = self.model.Part(
            name='Clamp', dimensionality=THREE_D, type=DISCRETE_RIGID_SURFACE)
        self.clampPart.BaseSolidExtrude(sketch=sketch, depth=4.0)
        RFPointID = self.clampPart.ReferencePoint(
            point=(coords[0], coords[1], 0.0)).id
        RFPoint = self.clampPart.referencePoints[RFPointID]
        RFPointRegion = self.clampPart.Set(referencePoints=(RFPoint,),
                                           name='Clamp Reference Point')
        # Create shell
        cells = self.clampPart.cells
        self.clampPart.RemoveCells(cellList=cells[0:1])
        # Assign inertia/mass to clamps
        self.clampPart.engineeringFeatures.PointMassInertia(
            name='Clamp Mass-Inertia', region=RFPointRegion, mass=1.0,
            alpha=0.0, composite=0.0)
        # Mesh clamp
        faces = self.clampPart.faces
        pickedRegions = (faces, )
        elemType1 = mesh.ElemType(elemCode=R3D4, elemLibrary=STANDARD)
        elemType2 = mesh.ElemType(elemCode=R3D3, elemLibrary=STANDARD)
        self.clampPart.setElementType(regions=pickedRegions,
                                      elemTypes=(elemType1, elemType2))
        self.clampPart.seedPart(size=1.0, deviationFactor=0.1,
                                minSizeFactor=0.1)
        self.clampPart.generateMesh()

    def createBottomPlate(self):
        '''Create the plate on which specimens rest for DWT simulations'''
        sketch = self.model.ConstrainedSketch(name='Bottom Plate Sketch',
                                              sheetSize=200.0)
        sketch.rectangle(point1=(-75.0, -100.0), point2=(75.0, 100.0))
        sketch.rectangle(point1=(-37.5, -62.5), point2=(37.5, 62.5))
        self.bPlatePart = self.model.Part(
            name='Bottom Plate', dimensionality=THREE_D,
            type=DISCRETE_RIGID_SURFACE)
        self.bPlatePart.BaseShell(sketch=sketch)
        RFPointID = self.bPlatePart.ReferencePoint(point=(0.0, 0.0, 0.0)).id
        RFPoint = self.bPlatePart.referencePoints[RFPointID]
        self.bPlatePart.Set(referencePoints=(RFPoint,),
                            name='Bottom Plate Reference Point')
        # Do not need to assign intertia as plate is constrained in every DOF
        # Mesh bottom plate
        faces = self.bPlatePart.faces
        pickedRegions = (faces, )
        elemType1 = mesh.ElemType(elemCode=R3D4, elemLibrary=STANDARD)
        elemType2 = mesh.ElemType(elemCode=R3D3, elemLibrary=STANDARD)
        self.bPlatePart.setElementType(regions=pickedRegions,
                                       elemTypes=(elemType1, elemType2))
        self.bPlatePart.seedPart(size=2.5, deviationFactor=0.1,
                                 minSizeFactor=0.1)
        self.bPlatePart.generateMesh()

    def createTapePaths(self, specimenDimensions):
        # create grids using shapely
        tapePaths, partGrid = tp.laminateCreation(
            tapeAngles=self.t_angles, tapeWidths=self.t_widths,
            tapeSpacing=self.t_spacing, specimenSize=specimenDimensions,
            undulationWidth=self.u_width)
        return tapePaths, partGrid

    def partition2DPart(self, part):
        '''
        Function to partition a 2D shell part into a grid of polygons
        depending on the tape angles of the laminate.
        '''
        # define mirror point (used to mirror the Polygon boundaries)
        mirrorPoint = Point([(0.0, 0.0), (0.0, 0.0)])

        tapes = zip(self.t_angles, self.t_widths)
        partitionLines = []
        uw = self.u_width
        for (a, w) in tapes:
            if a == 90:
                offset = w
                maxOffset = 50.0 - w / 2.0
                numberOffsets = int(maxOffset / offset)
                offsetList = [offset * k for k
                              in range(0, numberOffsets + 5)]
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
                reflectedLines = [scale(line, xfact=-1, origin=mirrorPoint)
                                  for line in partitionLines]
                partitionLines = partitionLines + reflectedLines

            else:
                offset = w / cos(radians(a))
                maxOffset = (75.0 - (w / 2.0) * cos(radians(a))
                             + 50.0 * tan(radians(a)))
                numberOffsets = int(maxOffset / offset)
                offsetList = [offset * k for k in range(0, numberOffsets + 5)]
                for ofs in offsetList:
                    tapeLineCoords = [(-100.0, (w / 2.0) - uw + ofs),
                                      (100.0, (w / 2.0) - uw + ofs)]
                    resinLineCoords1 = ([(-100.0, (w / 2.0) + ofs),
                                         (100.0, (w / 2.0) + ofs)])
                    resinLineCoords2 = ([(-100.0, (w / 2.0) + ofs + uw),
                                         (100.0, (w / 2.0) + ofs + uw)])
                    rotPoint = Point([(0.0, ofs), (0.0, ofs)])
                    tapeLine = rotate(LineString(tapeLineCoords), a, rotPoint)
                    resinLine1 = rotate(LineString(resinLineCoords1), a,
                                        rotPoint)
                    resinLine2 = rotate(LineString(resinLineCoords2), a,
                                        rotPoint)
                    partitionLines.extend([tapeLine, resinLine1, resinLine2])
                reflectedLines = [rotate(line, 180.0, origin=mirrorPoint)
                                  for line in partitionLines]
                partitionLines = partitionLines + reflectedLines

        for line in partitionLines:
            lineCoords = [list(tup) + [0.0] for tup in list(line.coords)]
            finalCoord = [[lineCoords[-1][0], lineCoords[-1][1], 1.0]]
            lineCoords3D = lineCoords + finalCoord
            dPointIDs = []
            for crds in lineCoords3D:
                dPointIDs.append(part.DatumPointByCoordinate(coords=crds).id)
            dPoint1 = part.datums[dPointIDs[0]]
            dPoint2 = part.datums[dPointIDs[1]]
            dPoint3 = part.datums[dPointIDs[2]]
            dpbtp = part.DatumPlaneByThreePoints(point1=dPoint1,
                                                 point2=dPoint2,
                                                 point3=dPoint3).id
            partFaces = part.faces
            try:
                part.PartitionFaceByDatumPlane(
                    datumPlane=part.datums[dpbtp], faces=partFaces)
            except:
                continue

    def create2DTapePart(self, xInterior, yInterior, xExterior, yExterior,
                         partGrid, meshSize):
        '''
        This function is used to create a (single!) 2D shell Tape part using
        geometric info from Shapely.
        '''
        # Create sketch
        sketch = self.model.ConstrainedSketch(name='2D Part Sketch',
                                              sheetSize=200.0)
        sketch.rectangle(point1=(-xExterior, -yExterior),
                         point2=(xExterior, yExterior))
        sketch.rectangle(point1=(-xInterior, -yInterior),
                         point2=(xInterior, yInterior))
        part = self.model.Part(name='Shell Part', dimensionality=THREE_D,
                               type=DEFORMABLE_BODY)
        part.BaseShell(sketch=sketch)
        partFaces = part.faces

        self.partition2DPart(part)

        for objID, obj in partGrid[1].iteritems():
            objCentroid = obj.centroid.coords
            selectedFace = partFaces.findAt(((objCentroid[0][0],
                                            objCentroid[0][1], 0.0), ))
            region = regionToolset.Region(faces=selectedFace)
            layup = part.CompositeLayup(
                name='CompositeLayup-{}'.format(objID), description='',
                elementType=SHELL, offsetType=MIDDLE_SURFACE, symmetric=False,
                thicknessAssignment=FROM_SECTION)
            layup.Section(
                preIntegrate=OFF, integrationRule=SIMPSON, useDensity=OFF,
                thicknessType=UNIFORM, temperature=GRADIENT,
                poissonDefinition=DEFAULT)
            layup.ReferenceOrientation(
                orientationType=GLOBAL, angle=0.0, localCsys=None, axis=AXIS_3,
                fieldName='', additionalRotationType=ROTATION_NONE)
            for layerNumber in range(1, len(partGrid) + 1):
                layerObj = partGrid[layerNumber][objID]
                objType = layerObj.objectType
                objAngle = layerObj.angle[-1]
                if objType == 'Undulation':
                    interfaceAngle = abs(layerObj.angle[0] - layerObj.angle[1])
                    secAngle = (0, interfaceAngle)
                    material = 'Undulation-{}'.format(secAngle)
                else:
                    material = objType
                nameOfPly = 'Ply-{}'.format(layerNumber)
                layup.CompositePly(
                    suppressed=False, thicknessType=SPECIFY_THICKNESS,
                    material=material, additionalRotationType=ROTATION_NONE,
                    thickness=0.2, orientationType=SPECIFY_ORIENT, angle=0.0,
                    numIntPoints=3, orientationValue=objAngle, axis=AXIS_3,
                    additionalRotationField='', plyName=nameOfPly,
                    region=region)

        partFaces = part.faces
        pickedRegions = (partFaces, )
        elemType1 = mesh.ElemType(elemCode=S4R, hourglassControl=ENHANCED,
                                  secondOrderAccuracy=OFF,
                                  elemLibrary=EXPLICIT)
        elemType2 = mesh.ElemType(elemCode=S3R, elemLibrary=EXPLICIT)
        part.setElementType(regions=pickedRegions,
                            elemTypes=(elemType1, elemType2))
        part.seedPart(size=meshSize, deviationFactor=0.1, minSizeFactor=0.1)
        part.generateMesh()

        instanceName = 'Shell Part Instance'
        assembly = self.model.rootAssembly
        assembly.Instance(name=instanceName, part=part, dependent=ON)
        assembly.rotate(instanceList=(instanceName, ), angle=-90.0,
                        axisPoint=(0.0, 0.0, 0.0),
                        axisDirection=(1.0, 0.0, 0.0))

    def create3DTapePart(self, obj, objectID, sequence, meshSize, damage=True,
                         symmetric=False):
        '''
        This function is used to create Tape parts using geometric info from
        Shapely.
        '''
        uniqueID = '{}-{}'.format(obj.layer, objectID)
        x0, y0 = obj.exterior.xy
        abqCoords = zip(x0, y0)
        sketch = self.model.ConstrainedSketch(
            name='Sketch {}'.format(uniqueID), sheetSize=200.0)
        for ind in range(len(abqCoords) - 1):
            sketch.Line(point1=(abqCoords[ind][0], abqCoords[ind][1]),
                        point2=(abqCoords[ind + 1][0], abqCoords[ind + 1][1]))
        # extrude profile to create part
        part = self.model.Part(name='Part {}'.format(uniqueID),
                               dimensionality=THREE_D, type=DEFORMABLE_BODY)
        part.BaseSolidExtrude(sketch=sketch, depth=self.t_thickness)
        cells = part.cells  # all cells within part (only 1 cell)
        cellRegion = regionToolset.Region(cells=cells[:])
        if obj.objectType == 'Tape':
            if damage:
                secName = 'Tape-Damage-Section'
            else:
                secName = 'Tape-Elastic-Section'
            part.SectionAssignment(
                region=cellRegion, sectionName=secName,
                offset=0.0, offsetType=MIDDLE_SURFACE, offsetField='',
                thicknessAssignment=FROM_SECTION)
        elif obj.objectType == 'Undulation':
            if len(obj.angle) == 1:
                if damage:
                    secName = 'Undulation-Damage-Resin-Section'
                else:
                    secName = 'Undulation-Elastic-Resin-Section'
            else:
                interfaceAngle = abs(obj.angle[0] - obj.angle[1])
                secAngle = (0, interfaceAngle)
                if damage:
                    secName = 'Undulation-Damage-{}-Section'.format(secAngle)
                else:
                    secName = 'Undulation-Elastic-{}-Section'.format(secAngle)
            part.SectionAssignment(
                region=cellRegion, thicknessAssignment=FROM_SECTION,
                offset=0.0, offsetType=MIDDLE_SURFACE, offsetField='',
                sectionName=secName)
        else:
            part.SectionAssignment(
                region=cellRegion, sectionName='Resin-Elastic-Section',
                offsetType=MIDDLE_SURFACE, thicknessAssignment=FROM_SECTION,
                offset=0.0, offsetField='')
        part.MaterialOrientation(
            region=cellRegion, orientationType=SYSTEM, stackDirection=STACK_3,
            angle=obj.angle[-1], additionalRotationType=ROTATION_ANGLE,
            axis=AXIS_3, fieldName='', additionalRotationField='',
            localCsys=None)

        # mesh the part
        part.seedPart(size=meshSize)
        part.generateMesh()

        # create the part instance
        instanceName = 'Instance3D-{}-{}'.format(sequence, uniqueID)
        sequenceOffset = sequence * self.t_thickness * len(self.t_angles)
        zOffset = sequenceOffset + (obj.layer - 1) * self.t_thickness
        self.assembly.Instance(name=instanceName, part=part, dependent=ON)
        self.assembly.translate(instanceList=(instanceName, ),
                                vector=(0.0, 0.0, zOffset))
        self.assembly.rotate(instanceList=(instanceName, ), angle=-90.0,
                             axisPoint=(0.0, 0.0, 0.0),
                             axisDirection=(1.0, 0.0, 0.0))
        # create the symmetric part instance
        if symmetric:
            print 'True'
            # create the symmetric part instance
            self.assembly.Instance(name='Mirror-' + instanceName, part=part,
                                   dependent=ON)
            self.assembly.translate(
                instanceList=('Mirror-' + instanceName, ),
                vector=(0.0, 0.0,
                        self.l_thickness - zOffset - self.t_thickness))
            objCentroid = obj.centroid.coords
            rotPoint = (objCentroid[0][0], objCentroid[0][1],
                        (self.l_thickness - zOffset) - self.t_thickness / 2.0)
            self.assembly.rotate(
                instanceList=('Mirror-' + instanceName, ), angle=180.0,
                axisPoint=rotPoint, axisDirection=(1.0, 0.0, 0.0))           
            self.assembly.rotate(
                instanceList=('Mirror-' + instanceName, ), angle=-90.0,
                axisPoint=(0.0, 0.0, 0.0), axisDirection=(1.0, 0.0, 0.0))

    def create3DTapePartGroup(self, partGrid, meshSize, symmetric=False):
        for sequence in range(self.sequences):
            for layerNumber in range(1, len(partGrid) + 1):
                for objID, obj in partGrid[layerNumber].iteritems():
                    self.create3DTapePart(obj, objID, sequence, meshSize,
                                          symmetric)

    def merge3DTapes(self, partGrid, tapePaths, meshSize, symmetric=False):
        if symmetric:
            for mirrorBool in (True, False):
                for sequence in range(self.sequences):
                    for (pathNumber, tapePath) in enumerate(tapePaths):
                        self.mergeInstances(partGrid, pathNumber, tapePath,
                                            sequence, meshSize,
                                            mirror=mirrorBool)
        else:
            for sequence in range(self.sequences):
                    for (pathNumber, tapePath) in enumerate(tapePaths):
                        self.mergeInstances(partGrid, pathNumber, tapePath,
                                            sequence, meshSize, mirror=False)

    def mergeInstances(self, partGrid, pathNumber, tapePath, sequence,
                       meshSize, mirror):
        if mirror:
            searchString = 'Mirror-Instance3D-{}-'.format(sequence)
        else:
            searchString = 'Instance3D-{}-'.format(sequence)
        if len(tapePath) > 1:
            layer, gridID = tapePath[0].split('-')
            pathAngle = partGrid[int(layer)][int(gridID)].angle[-1]
            instanceList = [self.assembly.instances[searchString + str(inst)]
                            for inst in tapePath]
            if mirror:
                instName = 'Mirror-Pass3D-{}-{}'.format(sequence, pathNumber)
            else:
                instName = 'Pass-3D-{}-{}'.format(sequence, pathNumber)
            self.assembly.InstanceFromBooleanMerge(
                name=instName, instances=instanceList, domain=BOTH,
                keepIntersections=ON, originalInstances=DELETE, mergeNodes=ALL,
                nodeMergingTolerance=1e-06)
            passPart = self.model.parts[instName]
            passRegion = regionToolset.Region(cells=passPart.cells[:])
            passPart.MaterialOrientation(
                region=passRegion, orientationType=SYSTEM, axis=AXIS_3,
                localCsys=None, fieldName='', stackDirection=STACK_3,
                additionalRotationType=ROTATION_ANGLE, angle=pathAngle,
                additionalRotationField='')
            passPart.seedPart(size=meshSize)
            passPart.generateMesh()

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

    def createAssembly(self):
        # Create the part instances
        self.assembly.Instance(name='Impactor Instance',
                               part=self.impactorPart, dependent=ON)
        self.assembly.Instance(name='Bottom Plate Instance',
                               part=self.bPlatePart, dependent=ON)
        for n in range(1, 5):
            instanceName = 'Clamp Instance {}'.format(n)
            self.assembly.Instance(name=instanceName, part=self.clampPart,
                                   dependent=ON)
        self.assembly.rotate(
            instanceList=('Bottom Plate Instance', 'Clamp Instance 1',
                          'Clamp Instance 2', 'Clamp Instance 3',
                          'Clamp Instance 4'),
            axisPoint=(0.0, 0.0, 0.0), axisDirection=(1.0, 0.0, 0.0),
            angle=-90.0)
        self.assembly.translate(instanceList=('Impactor Instance', ),
                                vector=(0.0, 8.75 + self.l_thickness, 0.0))
        self.assembly.translate(instanceList=('Clamp Instance 2', ),
                                vector=(-87.5, 0.0, 0.0))
        self.assembly.translate(instanceList=('Clamp Instance 3', ),
                                vector=(-87.5, 0.0, 100.0))
        self.assembly.translate(instanceList=('Clamp Instance 4', ),
                                vector=(0.0, 0.0, 100.0))
        self.assembly.translate(
            instanceList=('Clamp Instance 1', 'Clamp Instance 2',
                          'Clamp Instance 3', 'Clamp Instance 4'),
            vector=(0.0, self.l_thickness, 0.0))

    def constrainBPlate(self):
        # Encastre bottom plate
        bPlateInstance = self.assembly.instances['Bottom Plate Instance']
        bPlateBCRegion = bPlateInstance.sets['Bottom Plate Reference Point']
        self.model.EncastreBC(name='Bottom Plate Encastre BC',
                              createStepName='Loading Step', localCsys=None,
                              region=bPlateBCRegion)

    def constrainClamps(self):
        # Fix clamps in all but y-direction
        instances = self.assembly.instances
        clampInstances = [instances['Clamp Instance {}'.format(a)]
                          for a in range(1, 5)]
        clampRFPoints = [b.referencePoints.values()[0] for b in clampInstances]
        clampRegion = self.assembly.Set(referencePoints=clampRFPoints,
                                        name='Clamp BC Region')
        self.model.DisplacementBC(
            name='Clamp Displacement BC', createStepName='Loading Step',
            region=clampRegion, u1=0.0, u2=UNSET, u3=0.0, ur1=0.0, ur2=0.0,
            ur3=0.0, amplitude=UNSET, fixed=OFF, distributionType=UNIFORM,
            fieldName='', localCsys=None)
        # Apply load to clamps
        self.model.SmoothStepAmplitude(
            name='Smoothing Amplitude', timeSpan=STEP,
            data=((0.0, 0.0), (1e-05, 1.0)))
        self.model.ConcentratedForce(
            name='Clamp Load BC', createStepName='Loading Step',
            region=clampRegion, cf2=-1.0, amplitude='Smoothing Amplitude',
            distributionType=UNIFORM, field='', localCsys=None)

    def constrainImpactor(self):
        # Fix impactor in all but y direction
        impactorInstance = self.assembly.instances['Impactor Instance']
        impactorBCRegion = impactorInstance.sets['Impactor Reference Point']
        self.model.DisplacementBC(
            name='Impactor Displacement BC', createStepName='Loading Step',
            region=impactorBCRegion, u1=0.0, u2=UNSET, u3=0.0, ur1=0.0,
            ur2=0.0, ur3=0.0, amplitude=UNSET, fixed=OFF, fieldName='',
            distributionType=UNIFORM, localCsys=None)
        # Apply velocity BC to impactor
        self.model.Velocity(
            name='ImpactorVelocity', region=impactorBCRegion, field='', 
            distributionType=MAGNITUDE, velocity1=0.0,
            velocity2=self.initialVelocity, velocity3=0.0, omega=0.0)

    def createInteractionProperty(self):
        alpha = 50
        E_33 = 11.6
        G_12 = 6.47
        K1 = (alpha * E_33) / self.t_thickness
        K2 = (alpha * G_12) / self.t_thickness
        N = 0.055
        S = 0.105
        GIc = 300.0e-6
        GIIc = 800.0e-6
        self.model.ContactProperty('Tangential')
        self.model.interactionProperties['Tangential'].TangentialBehavior(
            formulation=PENALTY, directionality=ISOTROPIC,
            slipRateDependency=OFF, pressureDependency=OFF,
            temperatureDependency=OFF, dependencies=0, table=((0.15, ), ),
            shearStressLimit=None, maximumElasticSlip=FRACTION, fraction=0.005,
            elasticSlipStiffness=None)
        self.model.interactionProperties['Tangential'].NormalBehavior(
            pressureOverclosure=HARD, allowSeparation=ON,
            constraintEnforcementMethod=DEFAULT)
        self.model.ContactProperty('Cohesive')
        self.model.interactionProperties['Cohesive'].CohesiveBehavior(
            defaultPenalties=OFF, table=((K1, K2, K2), ))
        self.model.interactionProperties['Cohesive'].Damage(
            criterion=QUAD_TRACTION, initTable=((N, S, S), ),
            useEvolution=ON, evolutionType=ENERGY, useMixedMode=ON,
            mixedModeType=BK, exponent=1.75, evolTable=((GIc, GIIc, GIIc), ))

    def createInteractions(self):
        # determine contacts
        nonSpecimenInstances = ['Impactor Instance', 'Bottom Plate Instance',
                                'Clamp Instance 1', 'Clamp Instance 2',
                                'Clamp Instance 3', 'Clamp Instance 4']
        specimenInstances = [inst for inst in self.assembly.instances
                             if inst.name not in nonSpecimenInstances]
        self.model.contactDetection(
            defaultType=CONTACT, interactionProperty='Cohesive',
            nameEachSurfaceFound=OFF, createUnionOfMasterSurfaces=ON,
            createUnionOfSlaveSurfaces=ON, searchDomain=specimenInstances,
            separationTolerance=0.0001)
        # create explicit general contact definition
        self.model.ContactExp(name='GC', createStepName='Initial')
        generalContact = self.model.interactions['GC']
        generalContact.includedPairs.setValuesInStep(
            stepName='Initial', useAllstar=ON)
        # define 'Tangential' as default behaviour
        generalContact.contactPropertyAssignments.appendInStep(
            stepName='Initial', assignments=((GLOBAL, SELF, 'Tangential'), ))
        # assign cohesive behaviour to contacting tape surfaces
        for interact in self.model.interactions.values():
            if interact.name == 'GC':
                continue
            else:
                masterName = interact.master[0]
                slaveName = interact.slave[0]
                masterSurf = self.assembly.surfaces[masterName]
                slaveSurf = self.assembly.surfaces[slaveName]
                generalContact.contactPropertyAssignments.appendInStep(
                    stepName='Initial', 
                    assignments=((masterSurf, slaveSurf, 'Cohesive'), ))
                del self.model.interactions['{}'.format(interact.name)]

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
    energy = 30.0  # impact energy to simulate

    # Material parameters
    tapeAngles = (0, 90)  # define angles of tapes
    tapeWidths = 25.0
    tapeSpacing = 1  # number of gaps between tapes in interlacing pattern
    tapeThickness = 0.18  # cured ply thickness e.g. 0.18125
    undulationRatio = 0.18  # ratio of undulation amplitude to length
    nLayers = 24

    # Mesh sizes
    fineMesh = 0.25
    mediumMesh = 1.0
    coarseMesh = 3.0

    # Define specimen dimensions
    # Fine mesh region
    xMin_f = -25.0
    xMax_f = -xMin_f
    yMin_f = -25.0
    yMax_f = -yMin_f
    fineRegion = Polygon([(xMin_f, yMin_f), (xMax_f, yMin_f), (xMax_f, yMax_f),
                          (xMin_f, yMax_f)])

    # Shell region
    xInterior = -25.0
    xExterior = -50.0
    yInterior = -25.0
    yExterior = 75.0
    exterior = [(xExterior, yExterior), (xExterior, -yExterior),
                (-xExterior, -yExterior), (-xExterior, yExterior),
                (xExterior, yExterior)]
    interior = [(xInterior, yInterior), (xInterior, -yInterior),
                (-xInterior, -yInterior), (-xInterior, yInterior),
                (xInterior, yInterior)][::-1]
    shellRegion = Polygon(exterior, [interior])

    # Create model
    mdl = ImpactModel('DWT_Test')
    mdl.setTestParameters(time, outputIntervals, energy)
    mdl.setSpecimenParameters(tapeAngles, tapeWidths, tapeSpacing,
                              tapeThickness, undulationRatio, nLayers)
    mdl.createMaterials()
    # mdl.createImpactorPart()
    # mdl.createClampParts()
    # mdl.createBottomPlate()

    paths_f, grid_f = mdl.createTapePaths(fineRegion)
    # paths_sh, grid_sh = mdl.createTapePaths(shellRegion)
    # mdl.create2DTapePart(xInterior, yInterior, xExterior, yExterior, grid_sh,
    #                      mediumMesh)
    mdl.create3DTapePartGroup(grid_f, mediumMesh, symmetric=True)
    mdl.merge3DTapes(grid_f, paths_f, mediumMesh, symmetric=True)
    # mdl.createExplicitStep()
    # mdl.createAssembly()
    # mdl.constrainBPlate()
    # mdl.constrainClamps()
    # mdl.constrainImpactor()
    # mdl.createInteractionProperty()
    # mdl.createInteractions()

    # mdl.createJob()
