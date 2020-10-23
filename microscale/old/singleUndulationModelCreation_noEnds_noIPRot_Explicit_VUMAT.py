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

Current model limitations:
- in the undulating region cohesive zones only defined between plies of 
  differing orientation, rather than between each ply (reasonable since
  delamination rarely occurs between plies of the same orientation)

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
    mdb.Model(
            name='A-{}R-{}H-{}'.format(','.join(str(h_u).split('.')),
            ','.join(str(ratio).split('.')),','.join(str(h_f).split('.'))),
            modelType=STANDARD_EXPLICIT)
    undulationModel = mdb.models['A-{}R-{}H-{}'.format(
            ','.join(str(h_u).split('.')), ','.join(str(ratio).split('.')),
            ','.join(str(h_f).split('.')))]
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
    E_11 = 140.40
    E_22 = E_33 = 11.61
    nu_12 = nu_13 = 0.289
    nu_23 = 0.298
    G_12 = 6.47
    G_23 = 4.38
    nu_21 = nu_12*(E_22/E_11)
    nu_31 = nu_13*(E_33/E_22)
    nu_32 = nu_23*(E_33/E_22)

    tapeMaterial = laminateModel.Material(name='Tape')
    tapeMaterial.Density(table=((1.59e-06, ), ))
    tapeMaterial.Depvar(deleteVar=20, n=24)
    tapeMaterial.UserMaterial(
            mechanicalConstants=(146.8, 11.6, 11.6, 0.289, 0.289, 0.298, 6.47, 
            6.47, 4.38, 2.61, 1.759, 0.055, 0.285, 0.105, 53.0, 0.1, 0.1,
            0.00075, 0.0025, 0.0035))

    # Because it is difficult to define an orientation which varies both
    # in-plane and out-of plane, the out-of-plane undulation is defined
    # geometrically, while the in-plane orientation is taken into account by
    # modifying the material properties (rotating them)

    # Compliance
    S11 = 1/E_11
    S22 = S33 = 1/E_22
    S12 = S13 = -nu_12/E_11
    S23 = -nu_32/E_33
    S21 = S12
    S31 = S13
    S32 = S23
    S66 = S55 = 1/G_12
    S44 = 1/G_23

    S_lamina = np.array(
            [(S11, S12, S13, 0.0, 0.0, 0.0),
            (S21, S22, S23, 0.0, 0.0, 0.0),
            (S31, S32, S33, 0.0, 0.0, 0.0),
            (0.0, 0.0, 0.0, S44, 0.0, 0.0),
            (0.0, 0.0, 0.0, 0.0, S55, 0.0),
            (0.0, 0.0, 0.0, 0.0, 0.0, S66)])

    C_lamina = np.linalg.inv(S_lamina)
    
    S_undulation = np.linalg.inv(rotateLaminaStiffness(C_lamina, theta1))
    S_warp = np.linalg.inv(rotateLaminaStiffness(C_lamina, theta2))

    matprops_warp = engineeringConstants(S_warp)
    matprops_undulation = engineeringConstants(S_undulation)

    # Create materials
    undulationModel.Material(name='CF Warp')
    undulationModel.materials['CF Warp'].Elastic(type=ENGINEERING_CONSTANTS, 
        table=((matprops_warp), ))
    undulationModel.materials['CF Warp'].Density(table=((1.59e-06, ), ))

    undulationModel.Material(name='CF Undulation')
    undulationModel.materials['CF Undulation'].Elastic(
            type=ENGINEERING_CONSTANTS, table=((matprops_undulation), ))

    undulationModel.materials['CF Undulation'].Density(table=((1.59e-06, ), ))

    # ------------------------------------------------------------------------
    # Create sections to assign to the undulation
    cf_warp_section = undulationModel.HomogeneousSolidSection(
            name='CF Warp Section', material='CF Warp')
    cf_undulation_section = undulationModel.HomogeneousSolidSection(
            name='CF Undulation Section', material='CF Undulation')

    all_undulation_cells = undulationPart.Set(
            cells=undulationCells[:], name='All Undulation Cells')
    all_warp_cells = warpPart.Set(cells=warpCells[:], name='All Warp Cells')

    undulationPart.SectionAssignment(
            region=all_undulation_cells, sectionName='CF Undulation Section',
            offset=0.0, offsetType=MIDDLE_SURFACE, offsetField='',
            thicknessAssignment=FROM_SECTION)
    warpPart.SectionAssignment(
            region=all_warp_cells, sectionName='CF Warp Section', offset=0.0,
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
            axis=AXIS_1, normalAxisDefinition=VECTOR, 
            normalAxisVector=(0.0, 0.0, 1.0), flipNormalDirection=False, 
            normalAxisDirection=AXIS_3, primaryAxisDefinition=EDGE, 
            primaryAxisRegion=primaryAxisRegion, primaryAxisDirection=AXIS_1, 
            flipPrimaryDirection=False, additionalRotationType=ROTATION_NONE, 
            angle=0.0, additionalRotationField='', stackDirection=STACK_3)

    #select the two warp regions and assign material orientation
    warp_cell_region = regionToolset.Region(cells=warpCells[:])
    warpPart.MaterialOrientation(
            region=warp_cell_region, orientationType=SYSTEM, axis=AXIS_2,
            localCsys=None, fieldName='', 
            additionalRotationType=ROTATION_ANGLE, additionalRotationField='',
            angle=0.0, stackDirection=STACK_3)

    # <<<<<<<<<<<<<<<<<<<<<<<<<<< CREATE ASSEMBLY >>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    # Create the part instances
    undulationAssembly = undulationModel.rootAssembly
    undulationInstance = undulationAssembly.Instance(
            name='Undulation Instance', part=undulationPart, dependent=ON)

    warpInstance = undulationAssembly.Instance(
            name='Warp Instance', part=warpPart, dependent=ON)

    # <<<<<<<<<<<<<<<<<<<<<<<<<<<<< CREATE STEP >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    # create implicit step 
    undulationModel.ExplicitDynamicsStep(name='Loading Step',
            previous='Initial', timePeriod=2.5)

    undulationModel.fieldOutputRequests['F-Output-1'].setValues(
        variables=('S', 'SVAVG', 'LE', 'U', 'V', 'A', 'RF', 'CSTRESS',
            'CSDMG'))

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
            region=left_region, u1=0.0, u2=0.0, u3=0.0, ur1=0.0, ur2=0.0,
            ur3=0.0, amplitude=UNSET, fixed=OFF, distributionType=UNIFORM,
            fieldName='', localCsys=None)

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
    undulationModel.ContactProperty('Cohesive') # use default stiffnesses
    undulationModel.interactionProperties['Cohesive'].CohesiveBehavior(
        defaultPenalties=OFF, table=((K, K, K), ))
    undulationModel.interactionProperties['Cohesive'].Damage(
        initTable=((0.055, 0.09, 0.09), ))

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
    
    # Specify wedge element types for sliver regions
    # Use C3D10 tet elements. These are quadratic formulation elements. 
    # Contact should work well provided the formulation is surface to surface
    # and not node to surface
    elemType1 = mesh.ElemType(elemCode=UNKNOWN_HEX, elemLibrary=EXPLICIT)
    elemType2 = mesh.ElemType(elemCode=UNKNOWN_WEDGE, elemLibrary=EXPLICIT)
    elemType3 = mesh.ElemType(elemCode=C3D10M, elemLibrary=EXPLICIT)

    # Sliver regions
    thinCells = warpCells.findAt(((tan_p1-0.001, 0.001, w/2.0), ),
        ((tan_p2+0.001, 2*h_u-0.001, w/2.0), ))

    warpPart.setElementType(
            regions=(thinCells, ),
            elemTypes=(elemType1, elemType2, elemType3))
    warpPart.setMeshControls(
            regions=thinCells, elemShape=TET, technique=FREE)

    # Seed parts. Warp part must have a finer mesh than the undulation part for
    # correct contact enforcement
    warpPart.seedPart(size=0.015, deviationFactor=0.1, minSizeFactor=0.1)
    warpPart.generateMesh()

    undulationPart.seedPart(size=0.05, deviationFactor=0.1, minSizeFactor=0.1)
    undulationPart.generateMesh()    

    # <<<<<<<<<<<<<<<<<<<<<<<<<<<<< CREATE JOB >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Apply mass scaling to the thin geometry regions
    ms_cells = warpInstance.cells.findAt(((tan_p1-0.001, 0.001, w/2.0), ),
        ((tan_p2+0.001, 2*h_u-0.001, w/2.0), ))
    ms_set = undulationAssembly.Set(cells=ms_cells, name='MS Cells')
    massScalingRegion=undulationAssembly.sets['MS Cells']
    undulationModel.steps['Loading Step'].setValues(massScaling=((
        SEMI_AUTOMATIC, massScalingRegion, AT_BEGINNING, 16.0, 0.0, None, 0, 0, 0.0, 
        0.0, 0, None), ))
    
    # Create job and write input

    mdb.Job(name='A-{}R-{}H-{}'.format(','.join(str(h_u).split('.')),
            ','.join(str(ratio).split('.')),','.join(str(h_f).split('.'))),
            model='A-{}R-{}H-{}'.format(','.join(str(h_u).split('.')),
            ','.join(str(ratio).split('.')),','.join(str(h_f).split('.'))),
        description='', type=ANALYSIS, atTime=None, waitMinutes=0, waitHours=0,
        queue=None, memory=90, memoryUnits=PERCENTAGE,
        getMemoryFromAnalysis=True, explicitPrecision=DOUBLE_PLUS_PACK,
        nodalOutputPrecision=FULL, echoPrint=OFF, modelPrint=OFF,
        contactPrint=OFF, historyPrint=OFF, userSubroutine='', scratch='',
        resultsFormat=ODB, multiprocessingMode=DEFAULT, numCpus=1,
        numDomains=1, numGPUs=0, activateLoadBalancing=False)

    # Save model

    # mdb.saveAs(pathName='A-{}R-{}H-{}'.format(','.join(str(h_u).split('.')),
    #         ','.join(str(ratio).split('.')),','.join(str(h_f).split('.'))))

if __name__ == '__main__':
    # createUndulation(h_u, h_f, ratio, theta1, theta2):
    createUndulation(0.2,0.2,0.2, 0.0, 90.0)

    L_u = L_t = 1.0
    h_f = h_u = 0.2
    w = 0.1
    undulationModel = mdb.models['A-0,2R-0,2H-0,2']
    undulationPart = undulationModel.parts['Undulation']
    undulationEdges = undulationPart.edges
