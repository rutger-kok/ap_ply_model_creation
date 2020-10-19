from abaqus import *
from abaqusConstants import *
import mesh
import regionToolset
import os
from sys import path
githubPath = 'C:\\Users\\rutge\\Documents\\GitHub'
path.append('C:\\Python27\\Lib\\site-packages')
path.append(githubPath + '\\interlaced_model_creation\\editing')
from shapely.geometry import Polygon
import sigc3 as sigc
import interlacedMaterials as mats
import tapePlacement3 as tp


class TensileModel():
    def __init__(self, specimenType, undulationRatio):
        ratioString = str(undulationRatio).replace('.', '_')
        self.modelName = '{}-{}'.format(specimenType, ratioString)
        mdb.Model(name=self.modelName, modelType=STANDARD_EXPLICIT)
        self.model = mdb.models[self.modelName]
        self.assembly = self.model.rootAssembly
        # set the work directory
        wd = 'C:\\Workspace\\3D_Mechprops\\{}'.format(self.modelName)
        if not os.path.exists(wd):
            os.makedirs(wd)
        os.chdir(wd)

    def setTestParameters(self, time, outputIntervals):
        '''Set drop weight tower test parameters'''
        self.time = time
        self.outputIntervals = outputIntervals

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

    def createTapePaths(self, specimenDimensions):
        # create grids using shapely
        tapePaths, partGrid = tp.laminateCreation(
            tapeAngles=self.t_angles, tapeWidths=self.t_widths,
            tapeSpacing=self.t_spacing, specimenSize=specimenDimensions,
            undulationWidth=self.u_width, xMax=xMax, yMax=yMax)
        return tapePaths, partGrid

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
        if symmetric:
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

    def create3DTapePartGroup(self, partGrid, meshSize):
        for sequence in range(self.sequences):
            for layerNumber in range(1, len(partGrid) + 1):
                for objID, obj in partGrid[layerNumber].iteritems():
                    self.create3DTapePart(obj, objID, sequence, meshSize)

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
            region=rfPoint1Region, v1=UNSET, v2=0.05, v3=UNSET, vr1=UNSET,
            vr2=UNSET, vr3=UNSET, amplitude='Smoothing Amplitude', localCsys=None,
            distributionType=UNIFORM, fieldName='')
        self.model.ZsymmBC(name='Symmetry', createStepName='Loading Step',
                        region=symFacesSet, localCsys=None)

    def createInteractionProperty(self):
        alpha = 50
        E_33 = 9.10
        K = (alpha * E_33) / self.t_thickness

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
            defaultPenalties=OFF, table=((K, K, K), ))
        self.model.interactionProperties['Cohesive'].Damage(
            criterion=QUAD_TRACTION, initTable=((0.055, 0.09, 0.09), ),
            useEvolution=ON, evolutionType=ENERGY, useMixedMode=ON,
            mixedModeType=POWER_LAW, exponent=1.0,
            evolTable=((3.52e-4, 3.52e-4, 3.52e-4), ))

    def createInteractions(self):
        # determine contacts

        specimenInstances = [inst for inst in self.assembly.instances]
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


if __name__ == '__main__':

    # Simulation parameters
    time = 5.0  # duration to simulate [ms]
    outputIntervals = 50  # requested field output intervals

    # Material parameters
    specimenType = 'A'
    tapeAngles = (0, 90)  # define angles of tapes
    tapeWidths = 25.0
    tapeSpacing = 1  # number of gaps between tapes in interlacing pattern
    tapeThickness = 0.18  # cured ply thickness e.g. 0.18125
    undulationRatio = 0.09  # ratio of undulation amplitude to length
    nLayers = 4  # symmetric 8 ply laminate

    # Mesh sizes
    fineMesh = 0.25
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

    # Create model
    mdl = TensileModel(specimenType, undulationRatio)
    mdl.setTestParameters(time, outputIntervals)
    mdl.setSpecimenParameters(tapeAngles, tapeWidths, tapeSpacing,
                              tapeThickness, undulationRatio, nLayers)
    mdl.createMaterials()
    paths_f, grid_f = mdl.createTapePaths(fineRegion)
    mdl.create3DTapePartGroup(grid_f, mediumMesh)
    mdl.merge3DTapes(grid_f, paths_f, mediumMesh)
    mdl.createExplicitStep()
    mdl.applyContraints(dimensions)
    # mdl.createInteractionProperty()
    # mdl.createInteractions()

    # mdl.createJob()
