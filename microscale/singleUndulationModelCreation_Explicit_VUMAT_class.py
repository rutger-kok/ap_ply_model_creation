"""
This script file creates microscale models of fiber undulations.
Inputs to the model are:
- undulation amplitude h_u (mm)
- tape thickness h_f (mm)
- undulation ratio: ratio (-) 
Where the undulation ratio is the amplitude over the length of the undulation

The script creates the undulation profile in a sketch, uses this to partition
a box, splits the partitioned box into multiple parts (otherwise there can be 
no sliding between surfaces). 

Strengths and stiffnesses are given in GPa

"""

from abaqus import *
from abaqusConstants import *
import regionToolset
import mesh
import itertools
import numpy as np
from math import pi
from sys import path
githubPath = 'C:\\GitHub'
path.append('C:\\Python27\\Lib\\site-packages')
path.append(githubPath + '\\abaqus_single_undulation_model')

import interlacedMaterials2 as mats

class SingleUndulationModel():
    def __init__(self, h_u, ratio, theta1, theta2, w):
        ratioString = str(ratio).replace('.', '_')
        self.modelName = 'A{}-R{}-1{}-2{}'.format(specimenType, ratioString)
        mdb.Model(name=self.modelName, modelType=STANDARD_EXPLICIT)
        self.model = mdb.models[self.modelName]
        self.assembly = self.model.rootAssembly
        # set the work directory
        wd = 'C:\\Workspace\\3D_Mechprops\\{}'.format(self.modelName)
        if not os.path.exists(wd):
            os.makedirs(wd)
        os.chdir(wd)

        self.h_u = h_u
        self.ratio = ratio
        self.theta1 = theta1
        self.theta2 = theta2

        # define initial values
        self.L_u = self.h_u / self.ratio  # undulation length
        self.L_t = 1.0  # additional refined mesh region
        self.w = w  # specimen width
        self.L_g = self.L_u + 2 * self.L_t  # gauge length

    def setTestParameters(self, time, outputIntervals, testSpeed):
        '''Set test parameters'''
        self.time = time
        self.outputIntervals = outputIntervals
        self.testSpeed = testSpeed

    def createUndulationPart(self):
        # Create the undulating fiber part
        # Sketch rectangle
        undulationProfileSketch = self.model.ConstrainedSketch(
            name='Undulation Profile', sheetSize=5)
        undulationProfileSketch.rectangle(
            point1=(-self.L_t, 0.0),
            point2=(self.L_u + self.L_t, self.h_f + self.h_u))

        # Create a 3D deformable part named 'Undulation' by extruding the sketch
        # the part is simply a rectangular box
        undulationPart = self.model.Part(
            name='Undulation', dimensionality=THREE_D, type=DEFORMABLE_BODY)
        undulationPart.BaseSolidExtrude(
            sketch=undulationProfileSketch, depth=self.w)
        self.uPart = undulationPart

    def partitionUndulationPart(self, n=10):
        # Sketch the undulation and partition the part along its profile
        tangent = 1.0/6.0

        # Define tangents to curve to avoid excessively thin elements
        # linear curves y = mx + b
        tan_p1 = self.L_u*tangent  # point at which to calculate tangent curve
        grad_p1 = ((pi*self.h_u)/(2*self.L_u))*np.sin((pi*tan_p1)/self.L_u) # gradient
        ytan_p1 = self.h_u*(np.sin((pi*tan_p1)/(2*self.L_u)))**2  # y(x) at tan_p1
        btan_p1 = ytan_p1-grad_p1*tan_p1
        tan_p1_int = (0 - btan_p1)/grad_p1  # x coord at which y = mx+b = 0.0

        tan_p2 = self.L_u*(1-tangent)  # point at which to calculate tangent curve
        grad_p2 = ((pi*self.h_u)/(2*self.L_u))*np.sin((pi*tan_p2)/self.L_u) # gradient
        ytan_p2 = self.h_u*(np.sin((pi*tan_p2)/(2*self.L_u)))**2  # y(x) at tan_p2
        btan_p2 = ytan_p2-grad_p2*tan_p2
        tan_p2_int = (self.h_u - btan_p2)/grad_p2  # x coord at which y = mx+b = self.h_u
                
        self.profileCoords = [(0.0,0.0),(tan_p1_int,0.0),(tan_p1,ytan_p1)]  
        y_sine = lambda f: self.h_u*(np.sin((pi*f)/(2*self.L_u)))**2
        for x in np.linspace(tan_p1,tan_p2,n):
            x_coord = x
            y_coord = y_sine(x)
            self.profileCoords.append((x_coord,y_coord))
        self.profileCoords.extend([(tan_p2,ytan_p2),(tan_p2_int,self.h_u),(self.L_u,self.h_u)])

        profileUpper = [(a,b+self.h_u) for (a,b) in self.profileCoords] 
        undulationProfileSketch = self.model.ConstrainedSketch(
            name='Undulation Profile', sheetSize=5)
        
        for k in range(len(self.profileCoords)-1):
            undulationProfileSketch.Line(
                point1=self.profileCoords[k], point2=self.profileCoords[k+1])
            undulationProfileSketch.Line(
                point1=profileUpper[k], point2=profileUpper[k+1])

        # connect the upper and lower profile sketches to close the geometry
        undulationProfileSketch.Line(point1=(-self.L_t, 0.0), point2=(-self.L_t, self.h_f))
        undulationProfileSketch.Line(point1=(self.L_u+self.L_t, self.h_u), point2=(self.L_u+self.L_t, self.h_f+self.h_u))
        undulationProfileSketch.Line(point1=(-self.L_t, 0), point2=(0, 0))
        undulationProfileSketch.Line(point1=(self.L_u, self.h_u), point2=(self.L_u+self.L_t, self.h_u))
        undulationProfileSketch.Line(point1=(-self.L_t, self.h_f), point2=(0, self.h_f))
        undulationProfileSketch.Line(point1=(self.L_u, self.h_u+self.h_f), point2=(self.L_u+self.L_t, self.h_u+self.h_f))

        # define the partition sketch face
        partition_sketch_face = self.uPart.faces.findAt((self.L_u/2.0, self.h_u, self.w), )
        # choose an edge to define the orientation of the sketch partition
        # "choose an edge vertically and to the right"
        partition_edge = undulationEdges.findAt((self.L_u+self.L_t, self.h_u, self.w), )
        # select the cell to be partitioned
        partition_cell = undulationCells.findAt((self.L_u/2.0, self.h_u, self.w/2.0), )

        # partition the part by loading the previous undulation profile sketch
        partition_definition = self.uPart.MakeSketchTransform(
                sketchPlane=partition_sketch_face, sketchUpEdge=partition_edge,
                sketchPlaneSide=SIDE1, origin=(0.0, 0.0, self.w))
        partition_sketch = self.model.ConstrainedSketch(
                name='Partition Sketch', sheetSize=4.54, gridSpacing=0.11,
                transform=partition_definition)
        self.uPart.projectReferencesOntoSketch(
                sketch=partition_sketch, filter=COPLANAR_EDGES)
        partition_sketch.retrieveSketch(
                sketch=self.model.sketches['Undulation Profile'])
        self.uPart.PartitionCellBySketch(
                sketchPlane=partition_sketch_face, sketchUpEdge=partition_edge, 
                cells=partition_cell, sketch=partition_sketch)

        # for some reason the above changes the face definition so need to create 
        # another face for effective partitioning
        partition_face = self.uPart.faces.findAt((self.L_u/2.0, self.h_u, self.w), )  
        # select edges defining the partition profile
        partition_boundary = ([undulationEdges[b] for b in 
                partition_face.getEdges()])
        # select edges defining the partition path (sweep path)
        partition_path = undulationEdges.findAt((self.L_u+self.L_t, self.h_u+self.h_f, self.w/2.0), )
        self.uPart.PartitionCellBySweepEdge(
                sweepPath=partition_path, cells=partition_cell,
                edges=partition_boundary)

        # Cut sweep to remove warp parts of undulation
        self.uPart.CutSweep(path=undulationEdges.findAt(((-self.L_t, 2*self.h_u, self.w/4.0), )),
                profile=self.uPart.faces.findAt(coordinates=(0.0, 1.5*self.h_u, 
                    self.w)), flipSweepDirection=ON)

        self.uPart.CutSweep(path=undulationEdges.findAt(((-self.L_t, 0.0, self.w/4.0), )),
            profile=self.uPart.faces.findAt(coordinates=(self.L_u, 0.5*self.h_u, 
                self.w)), flipSweepDirection=OFF)


    def createWarpPart(self):
    
        # Copy the part and then cutsweep to create two separate parts: one for the
        # undulation and one for the two sections of warp tape. This is necessary
        # because otherwise the undulation and warp surfaces can not move relative
        # to eachother

        # copy undulation part
        self.wPart = self.model.Part(
                name='Warp', objectToCopy=self.uPart, compressFeatureList=ON)

        # remove undulating region from warp part
        self.wPart.CutSweep(
            path=self.wPart.edges.findAt(((-self.L_t, self.h_f+self.h_u, self.w/2.0), )),
            profile=self.wPart.faces.findAt(
                coordinates=(self.L_u/2.0, (self.h_f+self.h_u)/2.0, self.w)),
            flipSweepDirection=ON)

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
        mats.tapeDamage(self.model, E11, E22, E33, nu12, nu13, nu23, G12, G13,
                        G23, density, Xt, Xc, Yt, Yc, Sl, alpha0, G1Plus,
                        G1Minus, G2Plus, G2Minus, G6)

    def assignSection(self, damage=True):
        # Create sections to assign to the undulation
        if damage:
            tapeSection = self.model.HomogeneousSolidSection(
                name='Tape-Section', material='Tape-Damage')
        else:
            tapeSection = self.model.HomogeneousSolidSection(
                name='Tape-Section', material='Tape-Elastic')

        self.uCells = self.uPart.Set(
            cells=self.uPart.cells[:], name='All Undulation Cells')
        self.wCells = self.wPart.Set(
            cells=self.wPart.cells[:], name='All Warp Cells')

        self.uPart.SectionAssignment(
            region=alself.L_undulation_cells, sectionName='Tape-Section',
            offset=0.0, offsetType=MIDDLE_SURFACE, offsetField='',
            thicknessAssignment=FROM_SECTION)
        self.wPart.SectionAssignment(
            region=self.wCells, sectionName='Tape-Section', offset=0.0,
            offsetType=MIDDLE_SURFACE, offsetField='',
            thicknessAssignment=FROM_SECTION)

    def assignMatOrientation(self):
        # select the undulating cell and assign to region
        uPartCellRegion = regionToolset.Region(cells=self.uPart.cells[:])
        # create a set of the edges defining the material orientation 
        orientEdgeCoords = ([(((
            self.profileCoords[k][0]+self.profileCoords[k+1][0])/2.0,
            (self.profileCoords[k][1]+self.profileCoords[k+1][1])/2.0, self.w), ) 
            for k in range(len(self.profileCoords)-1)])
        orientCoordsString = str(orientEdgeCoords)[1:-1]
        findAtString = 'self.uPart.edges.findAt({})'.format(orientCoordsString)
        # need to use eval here because of specific findAt input requirements
        orientationEdges = eval(findAtString) 
        primaryAxisRegion = self.uPart.Set(
            edges=orientationEdges, name='Primary Axis Orientation')
        self.uPart.MaterialOrientation(
            region=uPartCellRegion, orientationType=DISCRETE,
            axis=AXIS_2, normalAxisDefinition=VECTOR, 
            normalAxisVector=(0.0, 0.0, 1.0), flipNormalDirection=False, 
            normalAxisDirection=AXIS_3, primaryAxisDefinition=EDGE, 
            primaryAxisRegion=primaryAxisRegion, primaryAxisDirection=AXIS_1, 
            flipPrimaryDirection=False, additionalRotationType=ROTATION_ANGLE, 
            angle=self.theta1, additionalRotationField='', stackDirection=STACK_3)

        #select the two warp regions and assign material orientation
        wPartCellRegion = regionToolset.Region(cells=self.wPart.cells[:])
        self.wPart.MaterialOrientation(
            region=wPartCellRegion, axis=AXIS_2, orientationType=SYSTEM,
            localCsys=None, fieldName='', angle=self.theta2, 
            additionalRotationType=ROTATION_ANGLE, additionalRotationField='',
            stackDirection=STACK_3)

    def meshParts(self, meshU=0.25, meshW=0.20):
        elemType1 = mesh.ElemType(elemCode=C3D8R, elemLibrary=EXPLICIT, 
            kinematicSplit=AVERAGE_STRAIN, secondOrderAccuracy=OFF, 
            hourglassControl=ENHANCED, distortionControl=DEFAULT)
        elemType2 = mesh.ElemType(elemCode=C3D6, elemLibrary=EXPLICIT)
        elemType3 = mesh.ElemType(elemCode=C3D4, elemLibrary=EXPLICIT)  

        self.wPart.seedPart(size=meshW, deviationFactor=0.1, minSizeFactor=0.1)
        self.wPart.generateMesh()
        self.wPart.setElementType(regions=self.wCells,
            elemTypes=(elemType1, elemType2, elemType3))

        self.uPart.seedPart(size=meshU, deviationFactor=0.1, minSizeFactor=0.1)
        self.uPart.generateMesh()    
        self.uPart.setElementType(regions=self.uCells,
            elemTypes=(elemType1, elemType2, elemType3))

    def createAssembly(self):
        # Create the part instances
        self.assembly.Instance(name='Undulation Instance', part=self.uPart,
                               dependent=ON)
        self.assembly.Instance(name='Warp Instance', part=self.wPart,
                               dependent=ON)

    def createExplicitStep(self):
        # mass scaling applied if increment < 1e-6
        self.model.ExplicitDynamicsStep(
            name='Loading Step', previous='Initial', timePeriod=self.time,
            massScaling=((SEMI_AUTOMATIC, MODEL, THROUGHOUT_STEP, 0.0, 1e-06,
                          BELOW_MIN, 1, 0, 0.0, 0.0, 0, None), ))
        self.model.fieldOutputRequests['F-Output-1'].setValues(
            variables=('S', 'LE', 'U', 'V', 'A', 'RF', 'CSTRESS', 'CSDMG',
                       'SDV', 'STATUS'), numIntervals=self.outputIntervals)

    def applyConstraints(self):
        self.model.SmoothStepAmplitude(name='Smoothing Amplitude',
            timeSpan=STEP, data=((0.0, 0.0), (1e-05, 1.0)))

        # Create partitions to define cells to apply BCs to
        offsetList = (-self.L_t + 0.25, 0, self.L_u,
                      self.L_u + self.L_t - 0.25, tan_p1, tan_p2)
        for part in self.model.parts.values():
            for offset in offsetList:
                dplane_id = part.DatumPlaneByPrincipalPlane(
                    principalPlane=YZPLANE, offset=offset).id
                dplane = part.datums[dplane_id]
                part.PartitionCellByDatumPlane(
                        datumPlane=dplane, cells=part.cells[:])

        d = 0.001
        leftBBox = (-self.L_t - d, -d, -d, -self.L_t + 0.26,
                    2.0 * self.h_u + d, d + self.w)
        leftCells = [inst.cells.getByBoundingBox(leftBBox) for inst
                     in self.assembly.instances.values() 
                     if inst.cells.getByBoundingBox(leftBBox)]
        leftRegion = self.assembly.Set(
            cells=left_cells, name='Left Coupling Region')
        rightBBox = (self.L_u + self.L_t - 0.26, -d, -d,
                     self.L_u + self.L_t + d, 2.0 * self.h_u + d, d + self.w)
        rightCells = [inst.cells.getByBoundingBox(rightBBox) for inst
                       in self.assembly.instances.values() 
                       if inst.cells.getByBoundingBox(rightBBox)]                
        rightRegion = self.assembly.Set(
            cells=rightCells, name='Right Coupling Region')    

        centerBBox = (-d, -d, -d, self.L_u + d,
                      2.0 * self.h_u + d, d + self.w)
        centerCells = [inst.cells.getByBoundingBox(centerBBox) for inst
                        in self.assembly.instances.values() 
                        if inst.cells.getByBoundingBox(centerBBox)]                
        centerRegion = self.assembly.Set(
            cells=centerCells, name='Center Cells')   

        # Apply displacement BC to left side of undulation
        self.model.DisplacementBC(
            name='Left Surface BC', createStepName='Loading Step',
            region=leftRegion, u1=0.0, u2=0.0, u3=0.0, ur1=UNSET,
            ur2=UNSET, ur3=UNSET, amplitude=UNSET, fixed=OFF,
            distributionType=UNIFORM, fieldName='', localCsys=None)

        # Apply velocity BC to right surface
        self.model.VelocityBC(name='Right Surface BC',
            createStepName='Loading Step', region=rightRegion,
            v1=self.testSpeed, v2=UNSET, v3=UNSET, vr1=UNSET, vr2=UNSET,
            vr3=UNSET, amplitude='Smoothing Amplitude', localCsys=None,
            distributionType=UNIFORM, fieldName='')

        # Apply symmetry across XZ plane
        symBBox = (-self.L_t - d, -d, -d, self.L_u + self.L_t + d,
                   d, d + self.w)
        symFaces = [inst.faces.getByBoundingBox(symBBox) for inst
                     in self.assembly.instances.values() 
                     if inst.faces.getByBoundingBox(symBBox)]   
        symRegion = self.assembly.Set(
            faces=symFaces, name='Symmetry Faces')
        # Symmetric across the XZ plane
        self.model.YsymmBC(name='Symmetry', 
            createStepName='Initial', region=symRegion, localCsys=None)

    def createInteractionProperty(self):
        alpha = 50
        E_33 = 11.6
        G_12 = 6.47
        GIc = 300.0e-6
        GIIc = 800.0e-6
        Ne = 5

        # calculate cohesive zone properties according to Turon (2007)
        K1 = (alpha * E_33) / self.t_thickness
        K2 = (alpha * G_12) / self.t_thickness
        tau1 = ((9*pi*E_33*GIc)/(32*Ne*self.meshSize))**0.5
        tau2 = ((9*pi*E_33*GIIc)/(32*Ne*self.meshSize))**0.5
        
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
            criterion=QUAD_TRACTION, initTable=((tau1, tau2, tau2), ),
            useEvolution=ON, evolutionType=ENERGY, useMixedMode=ON,
            mixedModeType=BK, exponent=1.75, evolTable=((GIc, GIIc, GIIc), ))

    def createInteractions(self):
        # determine contacts
        self.model.contactDetection(
            defaultType=CONTACT, interactionProperty='Cohesive',
            nameEachSurfaceFound=OFF, createUnionOfMasterSurfaces=ON,
            createUnionOfSlaveSurfaces=ON, searchDomain=MODEL,
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
    hu = 0.18
    hf = 0.18
    rat = 0.09
    width = 0.01
    phi1 = 0.0
    phi2 = 90.0

    mdl = SingleUndulationModel(hu, rat, phi1, phi2, width)
    mdl.setTestParameters(time=2.5, outputintervals=50, testSpeed=1.0)
    mdl.createUndulationPart()
    mdl.partitionUndulationPart()
    mdl.createWarpPart()
    mdl.createMaterials()
    mdl.assignSection()
    mdl.assignMatOrientation()
    mdl.meshParts()
    mdl.createAssembly()
    mdl.createExplicitStep()
    mdl.createInteractionProperty()
    mdl.createInteractions()
    mdl.applyConstraints()
    mdl.createJob()
    mdl.saveModel()


