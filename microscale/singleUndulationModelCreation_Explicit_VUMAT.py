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
from caeModules import *
import regionToolset
import mesh
import numpy as np
import math


def createUndulation(h_u, h_f, ratio, theta1, theta2):
    '''
    '''

    # define initial values
    L_u = h_u / ratio  # undulation length
    L_t = 1.0  # additional refined mesh region
    w = 0.1  # specimen width

    # initialize viewport and name model
    session.viewports['Viewport: 1'].setValues(displayedObject=None)
    modelName = 'A{}R{}H{}'.format(
        ''.join(str(h_u).split('.')), ''.join(str(ratio).split('.')),
        ''.join(str(h_f).split('.')))

    mdb.Model(name=modelName, modelType=STANDARD_EXPLICIT)
    undulationModel = mdb.models[modelName]
    undulationAssembly = undulationModel.rootAssembly

    # <<<<<<<<<<<<<<<<<<<<<<<<<<<< CREATE PARTS >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    # Create the undulating fiber part
    # Sketch rectangle
    undulationProfileSketch = undulationModel.ConstrainedSketch(
        name='Undulation Profile', sheetSize=5)
    undulationProfileSketch.rectangle(
        point1=(-L_t, 0.0), point2=(L_u + L_t, h_f + h_u))

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
    tangent = 1.0 / 6.0

    # Define tangents to curve to avoid excessively thin elements
    # linear curves y = mx + b
    tan_p1 = L_u * tangent  # point at which to calculate tangent curve
    grad_p1 = ((math.pi * h_u) / (2 * L_u)) * np.sin((math.pi * tan_p1) / L_u)
    ytan_p1 = h_u * (np.sin((math.pi * tan_p1) / (2 * L_u)))**2  # y(x) at tan_p1
    btan_p1 = ytan_p1 - grad_p1 * tan_p1
    tan_p1_int = (0 - btan_p1) / grad_p1  # x coord at which y = mx+b = 0.0

    tan_p2 = L_u * (1 - tangent)  # point at which to calculate tangent curve
    grad_p2 = ((math.pi * h_u) / (2 * L_u)) * np.sin((math.pi * tan_p2) / L_u)  # gradient
    ytan_p2 = h_u * (np.sin((math.pi * tan_p2) / (2 * L_u)))**2  # y(x) at tan_p2
    btan_p2 = ytan_p2 - grad_p2 * tan_p2
    tan_p2_int = (h_u - btan_p2) / grad_p2  # x coord at which y = mx+b = h_u

    profile_coords = [(0.0, 0.0), (tan_p1_int, 0.0), (tan_p1, ytan_p1)]
    y_sine = lambda f: h_u * (np.sin((math.pi * f) / (2 * L_u)))**2
    for x in np.linspace(tan_p1, tan_p2, n):
        x_coord = x
        y_coord = y_sine(x)
        profile_coords.append((x_coord, y_coord))
    profile_coords.extend([(tan_p2, ytan_p2), (tan_p2_int, h_u), (L_u, h_u)])

    profile_upper = [(a, b + h_u) for (a, b) in profile_coords]
    undulationProfileSketch = undulationModel.ConstrainedSketch(
        name='Undulation Profile', sheetSize=5)

    for k in range(len(profile_coords) - 1):
        undulationProfileSketch.Line(
            point1=profile_coords[k], point2=profile_coords[k + 1])
        undulationProfileSketch.Line(
            point1=profile_upper[k], point2=profile_upper[k + 1])

    # connect the upper and lower profile sketches to close the geometry
    undulationProfileSketch.Line(point1=(-L_t, 0.0), point2=(-L_t, h_f))
    undulationProfileSketch.Line(point1=(L_u + L_t, h_u), point2=(L_u + L_t, h_f + h_u))
    undulationProfileSketch.Line(point1=(-L_t, 0), point2=(0, 0))
    undulationProfileSketch.Line(point1=(L_u, h_u), point2=(L_u + L_t, h_u))
    undulationProfileSketch.Line(point1=(-L_t, h_f), point2=(0, h_f))
    undulationProfileSketch.Line(point1=(L_u, h_u + h_f), point2=(L_u + L_t, h_u + h_f))

    # define the partition sketch face
    partition_sketch_face = undulationPart.faces.findAt((L_u / 2.0, h_u, w), )
    # choose an edge to define the orientation of the sketch partition
    # "choose an edge vertically and to the right"
    partition_edge = undulationEdges.findAt((L_u + L_t, h_u, w), )
    # select the cell to be partitioned
    partition_cell = undulationCells.findAt((L_u / 2.0, h_u, w / 2.0), )

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
    partition_face = undulationPart.faces.findAt((L_u / 2.0, h_u, w), )
    # select edges defining the partition profile
    partition_boundary = (
        [undulationEdges[b] for b in partition_face.getEdges()])
    # select edges defining the partition path (sweep path)
    partition_path = undulationEdges.findAt((L_u + L_t, h_u + h_f, w / 2.0), )
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

    undulationPart.CutSweep(
        path=undulationEdges.findAt(((-L_t, 2 * h_u, w / 4.0), )),
        profile=undulationPart.faces.findAt(coordinates=(0.0, 1.5 * h_u, w)),
        flipSweepDirection=ON)

    undulationPart.CutSweep(
        path=undulationEdges.findAt(((-L_t, 0.0, w / 4.0), )),
        profile=undulationPart.faces.findAt(coordinates=(L_u, 0.5 * h_u, w)),
        flipSweepDirection=OFF)

    # remove undulating region from warp part
    warpFaces = warpPart.faces
    warpPart.CutSweep(
        path=warpPart.edges.findAt(((-L_t, h_f + h_u, w / 2.0), )),
        profile=warpFaces.findAt(coordinates=(L_u / 2.0, (h_f + h_u) / 2.0, w)),
        flipSweepDirection=ON)

    warpCells = warpPart.cells

    # <<<<<<<<<<<<<<<<<<<<<<<<<<<<< MATERIALS >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    # VTC401 - T700 UD estimated from Kitselis, Traiforos and Manolakos 2016
    # and Toray T700 datasheet
    # NOTE: these probably need to be changed when we get better experimental
    # values for the material properties

    tapeMaterial = undulationModel.Material(name='Tape')
    tapeMaterial.Density(table=((1.59e-06, ), ))
    tapeMaterial.Depvar(deleteVar=27, n=27)
    tapeMaterial.UserMaterial(
        mechanicalConstants=(116.6, 7.231, 7.231, 0.339, 0.339, 0.374, 3.268,
                             3.268, 2.632, 2.180, 0.2180, 0.811, 0.131, 0.185,
                             0.122, 0.070, 53.0, 0.5, 0.08, 0.02, 0.1, 0.00038,
                             0.00162))

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
    undulating_orientation_cell = regionToolset.Region(
        cells=undulationCells[:])
    # create a set of the edges defining the material orientation
    orient_edge_coords = (
        [(((profile_coords[k][0] + profile_coords[k + 1][0]) / 2.0,
           (profile_coords[k][1] + profile_coords[k + 1][1]) / 2.0, w), )
         for k in range(len(profile_coords) - 1)])
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

    # select the two warp regions and assign material orientation
    warp_cell_region = regionToolset.Region(cells=warpCells[:])
    warpPart.MaterialOrientation(
        region=warp_cell_region, orientationType=SYSTEM, axis=AXIS_2,
        localCsys=None, fieldName='', additionalRotationType=ROTATION_ANGLE,
        additionalRotationField='', angle=theta2, stackDirection=STACK_3)

    # <<<<<<<<<<<<<<<<<<<<<<<<<<< CREATE ASSEMBLY >>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    # Create the part instances
    undulationAssembly = undulationModel.rootAssembly
    undulationInstance = undulationAssembly.Instance(
        name='Undulation Instance', part=undulationPart, dependent=ON)

    warpInstance = undulationAssembly.Instance(
        name='Warp Instance', part=warpPart, dependent=ON)

    # <<<<<<<<<<<<<<<<<<<<<<<<<<<<< CREATE STEP >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

    # create explicit step
    undulationModel.ExplicitDynamicsStep(
        name='Loading Step', previous='Initial', timePeriod=5.0,
        massScaling=((SEMI_AUTOMATIC, MODEL, THROUGHOUT_STEP, 0.0,
                      1e-06, BELOW_MIN, 1, 0, 0.0, 0.0, 0, None), ))

    undulationModel.fieldOutputRequests['F-Output-1'].setValues(
        variables=('S', 'LE', 'U', 'V', 'A', 'RF', 'CSTRESS',
                   'CSDMG', 'SDV', 'STATUS'), numIntervals=50)

    # <<<<<<<<<<<<<<<<<<<<<<<<<< APPLY CONSTRAINTS >>>>>>>>>>>>>>>>>>>>>>>>>>>>

    def get_faces_by_bounding_box(bounding_box):
        '''
        This method searches every instance in an assembly to determine which
        faces are contained within a bounding box.

        Args:
            bounding_box (list): A list of coodinates defining a bounding box
                in the form: (x_min, y_min, z_min, x_max, y_max, z_max)

        Returns:
            faces (list): A list of the faces within a bounding box.
        '''
        x1, y1, z1, x2, y2, z2 = bounding_box
        faces = [inst.faces.getByBoundingBox(x1, y1, z1, x2, y2, z2)
                 for inst in undulationAssembly.instances.values()
                 if inst.faces.getByBoundingBox(x1, y1, z1, x2, y2, z2)]
        return faces

    def get_cells_by_bounding_box(bounding_box):
        '''
        This method searches every instance in an assembly to determine which
        cells are contained within a bounding box.

        Args:
            bounding_box (list): A list of coodinates defining a bounding box
                in the form: (x_min, y_min, z_min, x_max, y_max, z_max)

        Returns:
            cells (list): A list of the cells within a bounding box.
        '''
        x1, y1, z1, x2, y2, z2 = bounding_box
        cells = [inst.cells.getByBoundingBox(x1, y1, z1, x2, y2, z2)
                 for inst in undulationAssembly.instances.values()
                 if inst.cells.getByBoundingBox(x1, y1, z1, x2, y2, z2)]
        return cells

    # Create partitions to define cells to apply BCs to
    for prt in undulationModel.parts.values():
        for off_set in (-L_t + 0.25, 0, L_u, L_u + L_t - 0.25, tan_p1, tan_p2):
            dplane_id = prt.DatumPlaneByPrincipalPlane(
                principalPlane=YZPLANE, offset=off_set).id
            dplane = prt.datums[dplane_id]
            prt.PartitionCellByDatumPlane(
                datumPlane=dplane, cells=prt.cells[:])

    left_cells_coords = [-L_t - 0.001, -0.001, -0.001, -L_t + 0.26,
                         2.0 * h_u + 0.001, 0.001 + w]
    left_cells = get_cells_by_bounding_box(left_cells_coords)
    left_region = undulationAssembly.Set(
        cells=left_cells, name='Left Coupling Region')
    right_cells_coords = [L_u + L_t - 0.26, -0.001, -0.001, L_u + L_t + 0.001,
                          2.0 * h_u + 0.001, 0.001 + w]
    right_cells = get_cells_by_bounding_box(right_cells_coords)
    right_region = undulationAssembly.Set(
        cells=right_cells, name='Right Coupling Region')

    center_cells_coords = [-0.001, -0.001, -0.001, L_u + 0.001,
                           2.0 * h_u + 0.001, 0.001 + w]
    center_cells = get_cells_by_bounding_box(center_cells_coords)
    center_region = undulationAssembly.Set(
        cells=center_cells, name='Center Cells')

    # Apply displacement BC to left side of undulation
    undulationModel.DisplacementBC(
        name='Left Surface BC', createStepName='Loading Step',
        region=left_region, u1=0.0, u2=0.0, u3=0.0, ur1=UNSET,
        ur2=UNSET, ur3=UNSET, amplitude=UNSET, fixed=OFF,
        distributionType=UNIFORM, fieldName='', localCsys=None)

    # Apply velocity BC to right surface
    undulationModel.SmoothStepAmplitude(
        name='Smoothing Amplitude', timeSpan=STEP,
        data=((0.0, 0.0), (1e-05, 1.0)))
    undulationModel.VelocityBC(
        name='Right Surface BC', createStepName='Loading Step',
        region=right_region, v1=0.066, v2=UNSET, v3=UNSET, vr1=UNSET,
        vr2=UNSET, vr3=UNSET, amplitude='Smoothing Amplitude',
        localCsys=None, distributionType=UNIFORM, fieldName='')

    # Apply symmetry across XZ plane
    sym_faces_coords = [-L_t - 0.001, -0.001, -0.001, L_u + L_t + 0.001,
                        0.001, 0.001 + w]
    sym_faces = get_faces_by_bounding_box(sym_faces_coords)

    sym_region = undulationAssembly.Set(
        faces=sym_faces, name='Symmetry Faces')
    # Symmetric across the XZ plane
    undulationModel.YsymmBC(
        name='Symmetry', createStepName='Initial', region=sym_region,
        localCsys=None)

    undulationModel.SmoothStepAmplitude(
        name='Smoothing Amplitude', timeSpan=STEP,
        data=((0.0, 0.0), (1e-05, 1.0)))

    # <<<<<<<<<<<<<<<<<<<<<<<<< CREATE INTERACTIONS >>>>>>>>>>>>>>>>>>>>>>>>>>>

    # Create cohesive zone interactions between faces of the assembly
    # Calculate CZM parameters according to Turon (2007)
    alpha = 50  # from Turon (2007)
    E_33 = 11.6
    G_12 = 6.47
    GIc = 300.0e-6
    GIIc = 800.0e-6
    Ne = 5
    mesh_size = 0.05

    # calculate cohesive zone properties according to Turon (2007)
    K1 = (alpha * E_33) / h_u
    K2 = (alpha * G_12) / h_u
    tau1 = ((9 * pi * E_33 * GIc) / (32 * Ne * mesh_size))**0.5
    tau2 = ((9 * pi * E_33 * GIIc) / (32 * Ne * mesh_size))**0.5
    # Create cohesive zone interactions between faces of the assembly
    # In Abaqus Explicit the contact is defined as general contact. The default
    # behaviour is "hard" for normal contact, and penalty friction in the
    # tangential direction. Individual property assignments are used to
    # define the cohesive surfaces

    # define tangential material behaviour when not cohesive (friction)
    # current coefficient of friction set to 0.15
    undulationModel.ContactProperty('Tangential')
    undulationModel.interactionProperties['Tangential'].TangentialBehavior(
        formulation=PENALTY, directionality=ISOTROPIC,
        slipRateDependency=OFF, pressureDependency=OFF,
        temperatureDependency=OFF, dependencies=0, table=((0.15, ), ),
        shearStressLimit=None, maximumElasticSlip=FRACTION, fraction=0.005,
        elasticSlipStiffness=None)
    undulationModel.interactionProperties['Tangential'].NormalBehavior(
        pressureOverclosure=HARD, allowSeparation=ON,
        constraintEnforcementMethod=DEFAULT)
    undulationModel.ContactProperty('Cohesive')
    undulationModel.interactionProperties['Cohesive'].CohesiveBehavior(
        defaultPenalties=OFF, table=((K1, K2, K2), ))
    undulationModel.interactionProperties['Cohesive'].Damage(
        criterion=QUAD_TRACTION, initTable=((tau1, tau2, tau2), ),
        useEvolution=ON, evolutionType=ENERGY, useMixedMode=ON,
        mixedModeType=BK, exponent=1.75, evolTable=((GIc, GIIc, GIIc), ))

    # determine contacts
    undulationModel.contactDetection(
        defaultType=CONTACT, interactionProperty='Cohesive',
        nameEachSurfaceFound=OFF, createUnionOfMasterSurfaces=ON,
        createUnionOfSlaveSurfaces=ON,
        searchDomain=['Undulation Instance', 'Warp Instance'])

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

    elemType1 = mesh.ElemType(
        elemCode=C3D8R, elemLibrary=EXPLICIT, kinematicSplit=AVERAGE_STRAIN,
        secondOrderAccuracy=OFF, hourglassControl=ENHANCED,
        distortionControl=DEFAULT)
    elemType2 = mesh.ElemType(elemCode=C3D6, elemLibrary=EXPLICIT)
    elemType3 = mesh.ElemType(elemCode=C3D4, elemLibrary=EXPLICIT)  

    warpPart.seedPart(
        size=(mesh_size - 0.005), deviationFactor=0.1, minSizeFactor=0.1)
    warpPart.generateMesh()
    warpPart.setElementType(
        regions=all_warp_cells, elemTypes=(elemType1, elemType2, elemType3))

    undulationPart.seedPart(
        size=mesh_size, deviationFactor=0.1, minSizeFactor=0.1)
    undulationPart.generateMesh()
    undulationPart.setElementType(
        regions=all_undulation_cells,
        elemTypes=(elemType1, elemType2, elemType3))

    # <<<<<<<<<<<<<<<<<<<<<<<<<<<<< CREATE JOB >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    # Create job and write input

    # Define path to the material subroutine file (use a raw string)
    vumatPath = r"C:\GitHub\composite_cdm\composite_cdm.for"

    mdb.Job(
        name=modelName, model=modelName, description='', type=ANALYSIS,
        atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90,
        memoryUnits=PERCENTAGE, explicitPrecision=DOUBLE_PLUS_PACK,
        nodalOutputPrecision=FULL, echoPrint=OFF, modelPrint=OFF,
        contactPrint=OFF, historyPrint=OFF, userSubroutine=vumatPath,
        scratch='', resultsFormat=ODB, parallelizationMethodExplicit=DOMAIN,
        numDomains=6, activateLoadBalancing=False,
        multiprocessingMode=DEFAULT, numCpus=6)


if __name__ == '__main__':
    # createUndulation(h_u, h_f, ratio, theta1, theta2):
    # theta1 = in-plane orientation of undulating tape
    # theta2 = in-plane orientation of 'warp' section
    createUndulation(0.207, 0.207, 0.085, 0.0, 90.0)
