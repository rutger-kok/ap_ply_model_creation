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
import math

def createUndulation(h_u, h_f, ratio, theta1, theta2):
    '''
    '''

    # define initial values
    L_u = h_u/ratio  # undulation length
    L_t = 1.0 # additional refined mesh region
    w = 0.1 # specimen width
    L_g = L_u+2*L_t  # gauge length

    # initialize viewport and name model
    session.viewports['Viewport: 1'].setValues(displayedObject=None)
    modelName = 'A{}R{}H{}'.format(''.join(str(h_u).split('.')),
            ''.join(str(ratio).split('.')),''.join(str(h_f).split('.')))

    mdb.Model(name=modelName, modelType=STANDARD_EXPLICIT)
    undulationModel = mdb.models[modelName]
    undulationAssembly = undulationModel.rootAssembly

    # <<<<<<<<<<<<<<<<<<<<<<<<<<<< CREATE PARTS >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    
    # Create the undulating fiber part
    # Sketch rectangle
    undulationProfileSketch = undulationModel.ConstrainedSketch(
            name='Undulation Profile', sheetSize=5)
    undulationProfileSketch.rectangle(point1=(-L_t, 0.0),
            point2=(L_u+L_t, h_f+h_u))

    # Create a 3D deformable part named 'Undulation' by extruding the sketch
    # the part is simply a rectangular box
    undulationPart = undulationModel.Part(
            name='Undulation', dimensionality=THREE_D, type=DEFORMABLE_BODY)
    undulationPart.BaseSolidExtrude(sketch=undulationProfileSketch, depth=w)
    undulationCells = undulationPart.cells
    undulationEdges = undulationPart.edges

    # ------------------------------------------------------------------------
    # Sketch the undulation and partition the part along its profile
    n = 10  # number of points at which the undulation is discretized
    tangent = 1.0/6.0

    # Define tangents to curve to avoid excessively thin elements
    # linear curves y = mx + b
    tan_p1 = L_u*tangent  # point at which to calculate tangent curve
    grad_p1 = ((math.pi*h_u)/(2*L_u))*np.sin((math.pi*tan_p1)/L_u) # gradient
    ytan_p1 = h_u*(np.sin((math.pi*tan_p1)/(2*L_u)))**2  # y(x) at tan_p1
    btan_p1 = ytan_p1-grad_p1*tan_p1
    tan_p1_int = (0 - btan_p1)/grad_p1  # x coord at which y = mx+b = 0.0

    tan_p2 = L_u*(1-tangent)  # point at which to calculate tangent curve
    grad_p2 = ((math.pi*h_u)/(2*L_u))*np.sin((math.pi*tan_p2)/L_u) # gradient
    ytan_p2 = h_u*(np.sin((math.pi*tan_p2)/(2*L_u)))**2  # y(x) at tan_p2
    btan_p2 = ytan_p2-grad_p2*tan_p2
    tan_p2_int = (h_u - btan_p2)/grad_p2  # x coord at which y = mx+b = h_u
            
    profile_coords = [(0.0,0.0),(tan_p1_int,0.0),(tan_p1,ytan_p1)]  
    y_sine = lambda f: h_u*(np.sin((math.pi*f)/(2*L_u)))**2
    for x in np.linspace(tan_p1,tan_p2,n):
        x_coord = x
        y_coord = y_sine(x)
        profile_coords.append((x_coord,y_coord))
    profile_coords.extend([(tan_p2,ytan_p2),(tan_p2_int,h_u),(L_u,h_u)])

    profile_upper = [(a,b+h_u) for (a,b) in profile_coords] 
    undulationProfileSketch = undulationModel.ConstrainedSketch(
            name='Undulation Profile', sheetSize=5)
    
    for k in range(len(profile_coords)-1):
        undulationProfileSketch.Line(
                point1=profile_coords[k], point2=profile_coords[k+1])
        undulationProfileSketch.Line(
                point1=profile_upper[k], point2=profile_upper[k+1])

    # connect the upper and lower profile sketches to close the geometry
    undulationProfileSketch.Line(point1=(-L_t, 0.0), point2=(-L_t, h_f))
    undulationProfileSketch.Line(point1=(L_u+L_t, h_u), point2=(L_u+L_t, h_f+h_u))
    undulationProfileSketch.Line(point1=(-L_t, 0), point2=(0, 0))
    undulationProfileSketch.Line(point1=(L_u, h_u), point2=(L_u+L_t, h_u))
    undulationProfileSketch.Line(point1=(-L_t, h_f), point2=(0, h_f))
    undulationProfileSketch.Line(point1=(L_u, h_u+h_f), point2=(L_u+L_t, h_u+h_f))

    # define the partition sketch face
    partition_sketch_face = undulationPart.faces.findAt((L_u/2.0, h_u, w), )
    # choose an edge to define the orientation of the sketch partition
    # "choose an edge vertically and to the right"
    partition_edge = undulationEdges.findAt((L_u+L_t, h_u, w), )
    # select the cell to be partitioned
    partition_cell = undulationCells.findAt((L_u/2.0, h_u, w/2.0), )

    # partition the part by loading the previous undulation profile sketch
    partition_definition = undulationPart.MakeSketchTransform(
            sketchPlane=partition_sketch_face, sketchUpEdge=partition_edge,
            sketchPlaneSide=SIDE1, origin=(0.0, 0.0, w))
    partition_sketch = undulationModel.ConstrainedSketch(
            name='Partition Sketch', sheetSize=4.54, gridSpacing=0.11,
            transform=partition_definition)
    undulationPart.projectReferencesOntoSketch(
            sketch=partition_sketch, filter=COPLANAR_EDGES)
    partition_sketch.retrieveSketch(
            sketch=undulationModel.sketches['Undulation Profile'])
    undulationPart.PartitionCellBySketch(
            sketchPlane=partition_sketch_face, sketchUpEdge=partition_edge, 
            cells=partition_cell, sketch=partition_sketch)

    # for some reason the above changes the face definition so need to create 
    # another face for effective partitioning
    partition_face = undulationPart.faces.findAt((L_u/2.0, h_u, w), )  
    # select edges defining the partition profile
    partition_boundary = ([undulationEdges[b] for b in 
            partition_face.getEdges()])
    # select edges defining the partition path (sweep path)
    partition_path = undulationEdges.findAt((L_u+L_t, h_u+h_f, w/2.0), )
    undulationPart.PartitionCellBySweepEdge(
            sweepPath=partition_path, cells=partition_cell,
            edges=partition_boundary)

    # ------------------------------------------------------------------------
    # Copy the part and then cutsweep to create two separate parts: one for the
    # undulation and one for the two sections of warp tape. This is necessary
    # because otherwise the undulation and warp surfaces can not move relative
    # to eachother

    # copy undulation part
    warpPart = undulationModel.Part(
            name='Warp', objectToCopy=undulationPart, compressFeatureList=ON)

    undulationPart.CutSweep(path=undulationEdges.findAt(((-L_t, 2*h_u, w/4.0), )),
            profile=undulationPart.faces.findAt(coordinates=(0.0, 1.5*h_u, 
                w)), flipSweepDirection=ON)

    undulationPart.CutSweep(path=undulationEdges.findAt(((-L_t, 0.0, w/4.0), )),
        profile=undulationPart.faces.findAt(coordinates=(L_u, 0.5*h_u, 
            w)), flipSweepDirection=OFF)

    # # remove undulating region from warp part
    warpPart.CutSweep(
            path=warpPart.edges.findAt(((-L_t, h_f+h_u, w/2.0), )),
            profile=warpPart.faces.findAt(
                coordinates=(L_u/2.0, (h_f+h_u)/2.0, w)),
            flipSweepDirection=ON)

    warpCells = warpPart.cells
    warpEdges = warpPart.edges

    # <<<<<<<<<<<<<<<<<<<<<<<<<<<<< MATERIALS >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    # VTC401 - T700 UD estimated from Kitselis, Traiforos and Manolakos 2016
    # and Toray T700 datasheet
    # NOTE: these probably need to be changed when we get better experimental
    # values for the material properties

    tapeMaterial = undulationModel.Material(name='Tape')
    tapeMaterial.Density(table=((1.59e-06, ), ))
    tapeMaterial.Depvar(deleteVar=20, n=24)
    tapeMaterial.UserMaterial(
            mechanicalConstants=(146.8, 11.6, 11.6, 0.289, 0.289, 0.298, 6.47, 
            6.47, 4.38, 2.61, 1.759, 0.055, 0.285, 0.105, 53.0, 0.1, 0.1, 0.00075, 
            0.0025, 0.0035))
  
    # ------------------------------------------------------------------------
    # Create sections to assign to the undulation
    tapeSection = undulationModel.HomogeneousSolidSection(
            name='Tape Section', material='Tape')

    all_undulation_cells = undulationPart.Set(
            cells=undulationCells[:], name='All Undulation Cells')
    all_warp_cells = warpPart.Set(cells=warpCells[:], name='All Warp Cells')

    undulationPart.SectionAssignment(
            region=all_undulation_cells, sectionName='Tape Section',
            offset=0.0, offsetType=MIDDLE_SURFACE, offsetField='',
            thicknessAssignment=FROM_SECTION)
    warpPart.SectionAssignment(
            region=all_warp_cells, sectionName='Tape Section', offset=0.0,
            offsetType=MIDDLE_SURFACE, offsetField='',
            thicknessAssignment=FROM_SECTION)

    # ------------------------------------------------------------------------
    # Assign material orientations

    # select the undulating cell and assign to region
    # undulating_cell = undulationPart.cells.findAt(((a/2.0, h_u/2.0, w/2.0), ))
    undulating_orientation_cell = regionToolset.Region(
            cells=undulationCells[:])
    # create a set of the edges defining the material orientation 
    orient_edge_coords = ([(((profile_coords[k][0]+profile_coords[k+1][0])/2.0,
            (profile_coords[k][1]+profile_coords[k+1][1])/2.0, w), ) 
            for k in range(len(profile_coords)-1)])
    orient_coords_string = str(orient_edge_coords)[1:-1]
    findAt_string = 'undulationEdges.findAt({})'.format(orient_coords_string)
    # need to use eval here because of specific findAt input requirements
    orientation_edges = eval(findAt_string) 
    primaryAxisRegion = undulationPart.Set(
            edges=orientation_edges, name='Primary Axis Orientation')
    undulationPart.MaterialOrientation(
            region=undulating_orientation_cell, orientationType=DISCRETE,
            axis=AXIS_2, normalAxisDefinition=VECTOR, 
            normalAxisVector=(0.0, 0.0, 1.0), flipNormalDirection=False, 
            normalAxisDirection=AXIS_3, primaryAxisDefinition=EDGE, 
            primaryAxisRegion=primaryAxisRegion, primaryAxisDirection=AXIS_1, 
            flipPrimaryDirection=False, additionalRotationType=ROTATION_ANGLE, 
            angle=theta1, additionalRotationField='', stackDirection=STACK_3)

    #select the two warp regions and assign material orientation
    warp_cell_region = regionToolset.Region(cells=warpCells[:])
    warpPart.MaterialOrientation(
            region=warp_cell_region, orientationType=SYSTEM, axis=AXIS_2,
            localCsys=None, fieldName='', 
            additionalRotationType=ROTATION_ANGLE, additionalRotationField='',
            angle=theta2, stackDirection=STACK_3)

    # <<<<<<<<<<<<<<<<<<<<<<<<<<< CREATE ASSEMBLY >>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    # Create the part instances
    undulationAssembly = undulationModel.rootAssembly
    undulationInstance = undulationAssembly.Instance(
            name='Undulation Instance', part=undulationPart, dependent=ON)

    warpInstance = undulationAssembly.Instance(
            name='Warp Instance', part=warpPart, dependent=ON)

    # <<<<<<<<<<<<<<<<<<<<<<<<<<<<< CREATE STEP >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    # create explicit step 
    undulationModel.ExplicitDynamicsStep(name='Loading Step',
            previous='Initial', timePeriod=2.5, 
            massScaling=((SEMI_AUTOMATIC, MODEL, THROUGHOUT_STEP, 0.0,
            1e-06, BELOW_MIN, 1, 0, 0.0, 0.0, 0, None), ))

    undulationModel.fieldOutputRequests['F-Output-1'].setValues(
        variables=('S', 'LE', 'U', 'V', 'A', 'RF', 'CSTRESS',
            'CSDMG','SDV','STATUS'), numIntervals=50)

    # <<<<<<<<<<<<<<<<<<<<<<<<<< APPLY CONSTRAINTS >>>>>>>>>>>>>>>>>>>>>>>>>>>>

    # Create partitions to define cells to apply BCs to
    for prt in undulationModel.parts.values():
        for off_set in (-L_t+0.25, 0, L_u, L_u+L_t-0.25,tan_p1,tan_p2):
            dplane_id = prt.DatumPlaneByPrincipalPlane(
                principalPlane=YZPLANE, offset=off_set).id
            dplane = prt.datums[dplane_id]
            prt.PartitionCellByDatumPlane(
                    datumPlane=dplane, cells=prt.cells[:])

    left_cells = [inst.cells.getByBoundingBox(
                    -L_t-0.001, -0.001, -0.001, -L_t+0.26,
                    2.0*h_u+0.001, 0.001+w) for inst
                in undulationAssembly.instances.values() 
                if inst.cells.getByBoundingBox(
                    -L_t-0.001, -0.001, -0.001, -L_t+0.26,
                    2.0*h_u+0.001, 0.001+w)]
    left_region = undulationAssembly.Set(
                cells=left_cells, name='Left Coupling Region')
    right_cells = [inst.cells.getByBoundingBox(
                    L_u+L_t-0.26, -0.001, -0.001, L_u+L_t+0.001,
                    2.0*h_u+0.001, 0.001+w) for inst
                in undulationAssembly.instances.values() 
                if inst.cells.getByBoundingBox(
                    L_u+L_t-0.26, -0.001, -0.001, L_u+L_t+0.001,
                    2.0*h_u+0.001, 0.001+w)]                
    right_region = undulationAssembly.Set(
                cells=right_cells, name='Right Coupling Region')    

    center_cells = [inst.cells.getByBoundingBox(
                    -0.001, -0.001, -0.001, L_u+0.001,
                    2.0*h_u+0.001, 0.001+w) for inst
                in undulationAssembly.instances.values() 
                if inst.cells.getByBoundingBox(
                    -0.001, -0.001, -0.001, L_u+0.001,
                    2.0*h_u+0.001, 0.001+w)]                
    center_region = undulationAssembly.Set(
                cells=center_cells, name='Center Cells')   

    # Apply displacement BC to left side of undulation
    undulationModel.DisplacementBC(
            name='Left Surface BC', createStepName='Loading Step',
            region=left_region, u1=0.0, u2=0.0, u3=0.0, ur1=UNSET,
            ur2=UNSET, ur3=UNSET, amplitude=UNSET, fixed=OFF,
            distributionType=UNIFORM, fieldName='', localCsys=None)

    # Apply velocity BC to right surface
    undulationModel.SmoothStepAmplitude(name='Smoothing Amplitude',
            timeSpan=STEP, data=((0.0, 0.0), (1e-05, 1.0)))
    undulationModel.VelocityBC(name='Right Surface BC',
            createStepName='Loading Step', region=right_region, v1=0.02,
            v2=UNSET, v3=UNSET, vr1=UNSET, vr2=UNSET, vr3=UNSET,
            amplitude='Smoothing Amplitude', localCsys=None,
            distributionType=UNIFORM, fieldName='')

    # Apply symmetry across XZ plane
    sym_faces = [inst.faces.getByBoundingBox(
                    -L_t-0.001, -0.001, -0.001, L_u+L_t+0.001,
                    0.001, 0.001+w) for inst
                in undulationAssembly.instances.values() 
                if inst.faces.getByBoundingBox(
                    -L_t-0.001, -0.001, -0.001, L_u+L_t+0.001,
                    0.001, 0.001+w)]   

    sym_region = undulationAssembly.Set(
            faces=sym_faces, name='Symmetry Faces')
    # Symmetric across the XZ plane
    undulationModel.YsymmBC(name='Symmetry', 
            createStepName='Initial', region=sym_region, localCsys=None)

    undulationModel.SmoothStepAmplitude(name='Smoothing Amplitude',
        timeSpan=STEP, data=((0.0, 0.0), (1e-05, 1.0)))
            
    # <<<<<<<<<<<<<<<<<<<<<<<<< CREATE INTERACTIONS >>>>>>>>>>>>>>>>>>>>>>>>>>>
    
    # Create cohesive zone interactions between faces of the assembly
    # Calculate CZM parameters according to Turon (2007)
    alpha = 50
    E_33 = tapeMaterial.userMaterial.mechanicalConstants[2]
    K = (alpha*E_33)/h_f

    # Create cohesive zone interactions between faces of the assembly
    # In Abaqus Explicit the contact is defined as general contact. The default
    # behaviour is "hard" for normal contact, and penalty friction in the
    # tangential direction. Individual property assignments are used to 
    # define the cohesive surfaces

    # define tangential material behaviour when not cohesive (friction)
    # current coefficient of friction set to 0.15
    undulationModel.ContactProperty('Tangential')
    undulationModel.interactionProperties['Tangential'].TangentialBehavior(
        formulation=PENALTY, directionality=ISOTROPIC, slipRateDependency=OFF, 
        pressureDependency=OFF, temperatureDependency=OFF, dependencies=0, 
        table=((0.15, ), ), shearStressLimit=None, maximumElasticSlip=FRACTION, 
        fraction=0.005, elasticSlipStiffness=None)
    undulationModel.interactionProperties['Tangential'].NormalBehavior(
        pressureOverclosure=HARD, allowSeparation=ON, 
        constraintEnforcementMethod=DEFAULT)
    undulationModel.ContactProperty('Cohesive')
    undulationModel.interactionProperties['Cohesive'].CohesiveBehavior(
        defaultPenalties=OFF, table=((K, K, K), ))
    undulationModel.interactionProperties['Cohesive'].Damage(
        criterion=QUAD_TRACTION, initTable=((0.055, 0.09, 0.09), ), 
        useEvolution=ON, evolutionType=ENERGY, useMixedMode=ON, 
        mixedModeType=POWER_LAW, exponent=1.0,
        evolTable=((3.52e-4, 3.52e-4, 3.52e-4), ))

    # determine contacts
    undulationModel.contactDetection(defaultType=CONTACT,
        interactionProperty='Cohesive', nameEachSurfaceFound=OFF,
        createUnionOfMasterSurfaces=ON, createUnionOfSlaveSurfaces=ON,
        searchDomain=['Undulation Instance','Warp Instance'])

    undulationModel.ContactExp(name='GC', createStepName='Initial')
    undulationModel.interactions['GC'].includedPairs.setValuesInStep(
        stepName='Initial', useAllstar=ON)
    undulationModel.interactions['GC'].contactPropertyAssignments.appendInStep(
            stepName='Initial', assignments=((GLOBAL, SELF, 'Tangential'), ))
    for interact in undulationModel.interactions.values():
        if interact.name == 'GC': 
            continue
        else:
            masterName = interact.master[0]
            slaveName = interact.slave[0]
            masterSurf = undulationAssembly.surfaces[masterName]
            slaveSurf = undulationAssembly.surfaces[slaveName]
            undulationModel.interactions['GC'].contactPropertyAssignments.appendInStep(
                stepName='Initial', 
                assignments=((masterSurf, slaveSurf, 'Cohesive'), ))
            del undulationModel.interactions['{}'.format(interact.name)]

    # <<<<<<<<<<<<<<<<<<<<<<<<<<<< CREATE MESH >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    
    # Seed parts. Warp part must have a finer mesh than the undulation part for
    # correct contact enforcement

    elemType1 = mesh.ElemType(elemCode=C3D8R, elemLibrary=EXPLICIT, 
        kinematicSplit=AVERAGE_STRAIN, secondOrderAccuracy=OFF, 
        hourglassControl=ENHANCED, distortionControl=DEFAULT)
    elemType2 = mesh.ElemType(elemCode=C3D6, elemLibrary=EXPLICIT)
    elemType3 = mesh.ElemType(elemCode=C3D4, elemLibrary=EXPLICIT)  

    warpPart.seedPart(size=0.01, deviationFactor=0.1, minSizeFactor=0.1)
    warpPart.generateMesh()
    warpPart.setElementType(regions=all_warp_cells,
        elemTypes=(elemType1, elemType2, elemType3))

    undulationPart.seedPart(size=0.015, deviationFactor=0.1, minSizeFactor=0.1)
    undulationPart.generateMesh()    
    undulationPart.setElementType(regions=all_undulation_cells,
        elemTypes=(elemType1, elemType2, elemType3))

    # <<<<<<<<<<<<<<<<<<<<<<<<<<<<< CREATE JOB >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>    
    # Create job and write input

    # Define path to the material subroutine file (use a raw string)
    # vumatPath = r'C:\\workspace\\catalanotti_FC_v2.for'

    # mdb.Job(name=modelName, model=modelName, description='', type=ANALYSIS,
    #         atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90,
    #         memoryUnits=PERCENTAGE, explicitPrecision=DOUBLE_PLUS_PACK,
    #         nodalOutputPrecision=FULL, echoPrint=OFF, modelPrint=OFF,
    #         contactPrint=OFF, historyPrint=OFF, userSubroutine=vumatPath,
    #         scratch='', resultsFormat=ODB, parallelizationMethodExplicit=DOMAIN,
    #         numDomains=4, activateLoadBalancing=False,
    #         multiprocessingMode=DEFAULT, numCpus=4)

    # Save model

    # mdb.saveAs(pathName=modelName)

if __name__ == '__main__':
    # createUndulation(h_u, h_f, ratio, theta1, theta2):
    # theta1 = in-plane orientation of undulating tape
    # theta2 = in-plane orientation of 'warp' section
    createUndulation(0.207,0.207,0.085,0.0,90.0)


