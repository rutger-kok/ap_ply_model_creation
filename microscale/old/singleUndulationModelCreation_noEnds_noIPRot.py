"""
This script file creates microscale models of fiber undulations.
Inputs to the model are:
- undulation amplitude h_u (mm)
- tape thickness h_f (mm)
- undulation ratio: ratio (-) 
Where the undulation ratio is the amplitude over the length of the undulation
Strengths and stiffnesses are given in GPa.

The script creates the undulation profile in a sketch, uses this to partition
a box, splits the partitioned box into multiple parts (otherwise there can be 
no sliding between surfaces).
Materials are defined using engineering constants. The material orientation is 
aletered to account for the in-plane rotation of the "warp" tapes and the 
out-of-plane orientation of the undulating regions. 
A displacement BC is applied at the left end, and a velocity BC at the right 
end. Both are applied to reference points coupled to surfaces. Symmetry is 
applied across the XZ plane on the bottom surfaces.
Cohesive surfaces are defined between the contacting surfaces of the tapes.


Limitations:
- in the undulating region cohesive zones are only defined between plies of 
  differing orientation, rather than between each ply. This is a reasonable 
  assumeption since delamination rarely occurs between plies of the same
  orientation.
"""

from abaqus import *
from abaqusConstants import *
import regionToolset
import mesh
import itertools
import numpy as np
import math

def createUndulation(h_u, h_f, ratio):
    '''
    This function creates the undulation models.
    Inputs:
    - undulation amplitude h_u (mm)
    - tape thickness h_f (mm)
    - undulation ratio: ratio (-) 
    Outputs:
    - Abaqus CAE model of single tape undulation
    '''

    # define initial values
    L_u = h_u/ratio  # undulation length
    L_t = 1.0 # additional refined mesh region
    w = 0.1 # specimen width
    L_e = 0.1 # unused in this script, you can ignore this

    # initialize viewport and name model
    session.viewports['Viewport: 1'].setValues(displayedObject=None)
    # the names of the files are generated based on the parameters 
    # defining the undulation geometry. However, you cant use periods 
    # in Abaqus model names so here they are all replaced with commas.
    mdb.Model(
            name='A-{}R-{}H-{}'.format(','.join(str(h_u).split('.')),
            ','.join(str(ratio).split('.')),','.join(str(h_f).split('.'))),
            modelType=STANDARD_EXPLICIT)
    undulationModel = mdb.models['A-{}R-{}H-{}'.format(
            ','.join(str(h_u).split('.')), ','.join(str(ratio).split('.')),
            ','.join(str(h_f).split('.')))]
    undulationAssembly = undulationModel.rootAssembly

    # ------------------------------------------------------------------------
    # Create the undulating fiber part

    # Sketch rectangle
    undulationProfileSketch = undulationModel.ConstrainedSketch(
            name='Undulation Profile', sheetSize=5)
    undulationProfileSketch.rectangle(point1=(-L_t, 0.0), point2=(L_u+L_t, h_f+h_u))

    # Create a 3D deformable part named 'Undulation' by extruding the sketch
    # the part is simply a rectangular box
    undulationPart = undulationModel.Part(
            name='Undulation', dimensionality=THREE_D, type=DEFORMABLE_BODY)
    undulationPart.BaseSolidExtrude(sketch=undulationProfileSketch, depth=w)
    undulationCells = undulationPart.cells
    undulationEdges = undulationPart.edges

    # ------------------------------------------------------------------------
    # Sketch the undulation and partition the part along its profile
    n = 50  # number of points at which the undulation is discretized
    x = np.linspace(0,L_u,n)
    
    # The shape of the undulation is defined by a sine function. This results
    # in sliver geometries. To alleviate the issue (to an extent) I have 
    # defined a tangent to the sine curve at a point (L_u*tangent). This
    # reduces how thin the sliver geometries become. You can play around with
    # the tangent variable to ease meshing if necessary.

    # linear curves y = mx + b
    tangent = 2.0/6.0

    tan_p1 = L_u*tangent  # point at which to calculate tangent curve
    grad_p1 = ((math.pi*h_u)/(2*L_u))*np.sin((math.pi*tan_p1)/L_u) # gradient
    ytan_p1 = h_u*(np.sin((math.pi*tan_p1)/(2*L_u)))**2  # y(x) at tan_p1
    btan_p1 = ytan_p1-grad_p1*tan_p1
    tan_p1_int = (0 - btan_p1)/grad_p1  # x coord at which y = mx+b = h_u

    tan_p2 = L_u*(1-tangent)  # point at which to calculate tangent curve
    grad_p2 = ((math.pi*h_u)/(2*L_u))*np.sin((math.pi*tan_p2)/L_u) # gradient
    ytan_p2 = h_u*(np.sin((math.pi*tan_p2)/(2*L_u)))**2  # y(x) at tan_p1
    btan_p2 = ytan_p2-grad_p2*tan_p2
    tan_p2_int = (h_u - btan_p2)/grad_p2  # x coord at which y = mx+b = 0
    
    # define cross section profile of undulation
    h_u_prof = np.piecewise(
            x, [(x>=0.0)*(x<tan_p1_int), (x>=tan_p1_int)*(x<=tan_p1),
            (x>tan_p1)*(x<=tan_p2), (x>=tan_p2)*(x<=tan_p2_int),
            (x>tan_p2_int)*(x<=L_u)], 
            [lambda x: 0.0,
            lambda x: grad_p1*x + btan_p1,
            lambda x: h_u*(np.sin((math.pi*x)/(2*L_u)))**2,
            lambda x: grad_p2*x + btan_p2,
            lambda x: h_u])

    undulationProfileSketch = undulationModel.ConstrainedSketch(
            name='Undulation Profile', sheetSize=5)

    # reduce number of edges by merging all y = 0 and y=h_u edges
    bottom_last = np.where(h_u_prof == 0)[0][-1]
    top_first = np.where(h_u_prof == h_u)[0][0]
    mod_x = np.concatenate(([x[0]],x[bottom_last:top_first+1],[x[-1]]))
    mod_h_u_prof = np.concatenate(
            ([h_u_prof[0]],h_u_prof[bottom_last:top_first+1],[h_u_prof[-1]]))
    
    mod_h_u_prof[mod_h_u_prof < 0.001] = 0.0  # eliminate thin sections
    mod_h_u_prof[(h_u-mod_h_u_prof) < 0.001] = h_u
    h_f_prof = mod_h_u_prof+h_f  # curve defining the top of the undulation

    for k in range(len(mod_x)-1):
        undulationProfileSketch.Line(
                point1=(mod_x[k], mod_h_u_prof[k]),
                point2=(mod_x[k+1], mod_h_u_prof[k+1]))
        undulationProfileSketch.Line(
                point1=(mod_x[k], h_f_prof[k]),
                point2=(mod_x[k+1], h_f_prof[k+1]))

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

    # delete the "warp" parts
    undulationPart.CutSweep(
            path=undulationEdges.findAt(((-L_t, 0.0, w/2.0), )),
            profile=undulationPart.faces.findAt(
                coordinates=(0.001, (h_u+h_f)-0.001, w)),
            flipSweepDirection=ON)  # top warp section
    undulationPart.CutSweep(
            path=undulationEdges.findAt(((-L_t, 0.0, w/2.0), )),
            profile=undulationPart.faces.findAt(
                coordinates=(L_u-0.001, 0.001, w)),
            flipSweepDirection=ON)  # bottom warp section

    # remove undulating region from warp part
    warpPart.CutSweep(
            path=warpPart.edges.findAt(((-L_t, h_f+h_u, w/2.0), )),
            profile=warpPart.faces.findAt(
                coordinates=(L_u/2.0, (h_f+h_u)/2.0, w)),
            flipSweepDirection=ON)

    warpCells = warpPart.cells
    warpEdges = warpPart.edges

    # ------------------------------------------------------------------------
    # Because it is difficult to define an orientation which varies both
    # in-plane and out-of plane, the out-of-plane undulation is defined
    # geometrically, while the in-plane orientation is taken into account by
    # modifying the material properties (rotating them)

    # VTC401 - T700 UD estimated from Kitselis, Traiforos and Manolakos 2016
    # and Toray T700 datasheet
    # NOTE: these probably need to be changed when we get better experimental
    # values for the material properties
    E_11 = 140.40
    E_22 = E_33 = 11.61
    nu_12 = nu_13 = 0.289
    nu_23 = 0.298
    G_12 = G_13 = 6.47
    G_23 = 4.38
    nu_21 = nu_12*(E_22/E_11)
    nu_31 = nu_13*(E_33/E_22)
    nu_32 = nu_23*(E_33/E_22)

    matProps = (E_11,E_22,E_33,nu_12,nu_13,nu_23,G_12,G_13,G_23)

    # Create materials
    undulationModel.Material(name='CF Tape')
    undulationModel.materials['CF Tape'].Elastic(type=ENGINEERING_CONSTANTS, 
        table=((matProps), ))
    undulationModel.materials['CF Tape'].Density(table=((1.59e-06, ), ))

    # ------------------------------------------------------------------------
    # Create sections to assign to the undulation
    cf_section = undulationModel.HomogeneousSolidSection(
            name='CF Section', material='CF Tape')

    all_undulation_cells = undulationPart.Set(
            cells=undulationCells[:], name='All Undulation Cells')
    all_warp_cells = warpPart.Set(cells=warpCells[:], name='All Warp Cells')

    undulationPart.SectionAssignment(
            region=all_undulation_cells, sectionName='CF Section',
            offset=0.0, offsetType=MIDDLE_SURFACE, offsetField='',
            thicknessAssignment=FROM_SECTION)
    warpPart.SectionAssignment(
            region=all_warp_cells, sectionName='CF Section', offset=0.0,
            offsetType=MIDDLE_SURFACE, offsetField='',
            thicknessAssignment=FROM_SECTION)

    # ------------------------------------------------------------------------
    # Assign material orientations

    # select the undulating cell and assign to region
    # undulating_cell = undulationPart.cells.findAt(((a/2.0, h_u/2.0, w/2.0), ))
    undulating_orientation_cell = regionToolset.Region(
            cells=undulationCells[:])
    # create a set of the edges defining the material orientation 
    orient_edge_coords = ([(((mod_x[k]+mod_x[k+1])/2.0,
            (mod_h_u_prof[k]+mod_h_u_prof[k+1])/2.0, w), ) 
            for k in range(len(mod_x)-1)])
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

    #select the two warp regions and assign material orientation (90 degrees)
    warp_cell_region = regionToolset.Region(cells=warpCells[:])
    warpPart.MaterialOrientation(
            region=warp_cell_region, orientationType=SYSTEM, axis=AXIS_2,
            localCsys=None, fieldName='',
            additionalRotationType=ROTATION_ANGLE, additionalRotationField='',
            angle=90.0, stackDirection=STACK_3)

    # ------------------------------------------------------------------------
    # Create the assembly

    # Create the part instance
    undulationAssembly = undulationModel.rootAssembly
    undulationInstance = undulationAssembly.Instance(
            name='Undulation Instance', part=undulationPart, dependent=ON)

    warpInstance = undulationAssembly.Instance(
            name='Warp Instance', part=warpPart, dependent=ON)

    # ------------------------------------------------------------------------
    # Create the step

    # create explicit step with duration 10s
    undulationModel.ExplicitDynamicsStep(name='Loading Step',
            previous='Initial', timePeriod=10.0)

    # modify field output request
    undulationModel.fieldOutputRequests['F-Output-1'].setValues(
            numIntervals=100)

    # ------------------------------------------------------------------------
    # Apply constraints

    # Create reference points on either side of the undulation
    # create coupling constraints between the sides of the undulation and the 
    # reference points
    rfPoint1_id = undulationAssembly.ReferencePoint(
            point=(-L_e-L_t-1, 0.0, w/2.0)).id
    rfPoint1 = undulationAssembly.referencePoints[rfPoint1_id]
    rfPoint1_region = undulationAssembly.Set(
            referencePoints=(rfPoint1,), name='Coupling Reference Point 1')

    # Create reference point to the right of the undulation
    rfPoint2_id = undulationAssembly.ReferencePoint(
            point=(L_u+L_e+L_t+1, 0.0, w/2.0)).id
    rfPoint2 = undulationAssembly.referencePoints[rfPoint2_id]
    rfPoint2_region = undulationAssembly.Set(
            referencePoints=(rfPoint2,), name='Coupling Reference Point 2')

    # Create coupling constraint between bottom surface and reference point
    left_coupling_face = [inst.faces.getByBoundingBox(
                    -L_t-0.001, -0.001, -0.001, -L_t+0.001,
                    (h_f+h_u)+0.001, 0.001+w) for inst
                in undulationAssembly.instances.values() 
                if inst.faces.getByBoundingBox(
                    -L_t-0.001, -0.001, -0.001, -L_t+0.001,
                    (h_f+h_u)+0.001, 0.001+w)]
    left_coupling_region = undulationAssembly.Surface(
            side1Faces=left_coupling_face,
            name='Left Coupling Surface')
    right_coupling_face = [inst.faces.getByBoundingBox(
                L_u+L_t-0.001, -0.001, -0.001, L_u+L_t+0.001,
                (h_f+h_u)+0.001, 0.001+w) for inst
            in undulationAssembly.instances.values() 
            if inst.faces.getByBoundingBox(
                L_u+L_t-0.001, -0.001, -0.001, L_u+L_t+0.001,
                (h_f+h_u)+0.001, 0.001+w)]
    right_coupling_region = undulationAssembly.Surface(
            side1Faces=right_coupling_face,
            name='Right Coupling Surface')
        
    undulationModel.Coupling(
            name='Left Surface', controlPoint=rfPoint1_region,
            surface=left_coupling_region, influenceRadius=WHOLE_SURFACE,
            couplingType=KINEMATIC, localCsys=None, u1=ON, u2=ON, u3=ON,
            ur1=ON, ur2=ON, ur3=ON)

    undulationModel.Coupling(
            name='Right Surface', controlPoint=rfPoint2_region,
            surface=right_coupling_region, influenceRadius=WHOLE_SURFACE,
            couplingType=KINEMATIC, localCsys=None, u1=ON, u2=ON, u3=ON,
            ur1=ON, ur2=ON, ur3=ON)

    # ------------------------------------------------------------------------
    # Apply BCs
    # Constrain left surface
    undulationModel.DisplacementBC(
            name='Left Surface BC', createStepName='Loading Step',
            region=rfPoint1_region, u1=0.0, u2=0.0, u3=0.0, ur1=0.0, ur2=0.0,
            ur3=0.0, amplitude=UNSET, fixed=OFF, distributionType=UNIFORM,
            fieldName='', localCsys=None)

    # Apply velocity BC to right surface
    undulationModel.SmoothStepAmplitude(name='Smoothing Amplitude',
            timeSpan=STEP, data=((0.0, 0.0), (1e-05, 1.0)))
    undulationModel.VelocityBC(name='Right Surface BC',
            createStepName='Loading Step', region=rfPoint2_region, v1=0.02,
            v2=UNSET, v3=UNSET, vr1=UNSET, vr2=UNSET, vr3=UNSET,
            amplitude='Smoothing Amplitude', localCsys=None,
            distributionType=UNIFORM, fieldName='')

    # Apply symmetry across the XZ plane
    sym_face1 = warpInstance.faces.findAt(((L_u-0.001, 0.0, w/2.0), ))
    sym_face2 = undulationInstance.faces.findAt(((0.001, 0.0, w/2.0), ))
    sym_region = undulationAssembly.Set(
            faces=sym_face1+sym_face2, name='Symmetry Faces')
    
    undulationModel.YsymmBC(name='Symmetry', 
            createStepName='Initial', region=sym_region, localCsys=None)
            
    # ------------------------------------------------------------------------
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
    undulationModel.interactionProperties['Cohesive'].CohesiveBehavior()
    # undulationModel.interactionProperties['Cohesive'].Damage(initTable=((
    #     1.0, 2.0, 3.0), ), evolTable=())

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
                stepName='Initial', assignments=((masterSurf, slaveSurf, 'Cohesive'), ))
            del undulationModel.interactions['{}'.format(interact.name)]

    # ------------------------------------------------------------------------
    # Create a partition and a set of nodes to define the undulation region
    # this is used to determine the stiffness of the undulation region only
    leftside_warp_dplane_id = warpPart.DatumPlaneByPrincipalPlane(
            principalPlane=YZPLANE, offset=0.0).id
    leftside_warp_dplane = warpPart.datums[leftside_warp_dplane_id]
    rightside_warp_dplane_id = warpPart.DatumPlaneByPrincipalPlane(
            principalPlane=YZPLANE, offset=L_u).id
    rightside_warp_dplane = warpPart.datums[rightside_warp_dplane_id]
    warpPart.PartitionCellByDatumPlane(
            datumPlane=leftside_warp_dplane, cells=warpCells[:])
    warpPart.PartitionCellByDatumPlane(
            datumPlane=rightside_warp_dplane, cells=warpCells[:])

    leftside_undulation_dplane_id = undulationPart.DatumPlaneByPrincipalPlane(
            principalPlane=YZPLANE, offset=0.0).id
    leftside_undulation_dplane = (
            undulationPart.datums[leftside_undulation_dplane_id])
    rightside_undulation_dplane_id = undulationPart.DatumPlaneByPrincipalPlane(
            principalPlane=YZPLANE, offset=L_u).id
    rightside_undulation_dplane = (
            undulationPart.datums[rightside_undulation_dplane_id])
    undulationPart.PartitionCellByDatumPlane(
            datumPlane=leftside_undulation_dplane, cells=undulationCells[:])
    undulationPart.PartitionCellByDatumPlane(
            datumPlane=rightside_undulation_dplane, cells=undulationCells[:])

    leftside_warp_faces = warpInstance.faces.findAt(
            ((0.0, h_f+h_u/2.0, w/2.0), ))
    leftside_undulation_faces = undulationInstance.faces.findAt(
            ((0.0, h_f/2.0, w/2.0), ))
    rightside_warp_faces = warpInstance.faces.findAt(
            ((L_u, h_u/2.0, w/2.0), ))
    rightside_undulation_faces = undulationInstance.faces.findAt(
            ((L_u, h_u+h_f/2.0, w/2.0), ))
    undulationAssembly.Set(
            faces=leftside_undulation_faces+leftside_warp_faces,
            name='Left Side')
    undulationAssembly.Set(
            faces=rightside_undulation_faces+rightside_warp_faces,
            name='Right Side')

    # ------------------------------------------------------------------------
    # Create mesh
    # You can modify the size argument here to change the mesh fineness

    undulationPart.seedPart(size=0.05, deviationFactor=0.1, minSizeFactor=0.1)
    undulationPart.generateMesh()

    warpPart.seedPart(size=0.05, deviationFactor=0.1, minSizeFactor=0.1)
    warpPart.generateMesh()

    # ------------------------------------------------------------------------
    # Create job

    mdb.Job(name='A-{}R-{}H-{}'.format(','.join(str(h_u).split('.')),
            ','.join(str(ratio).split('.')),','.join(str(h_f).split('.'))),
            model='A-{}R-{}H-{}'.format(','.join(str(h_u).split('.')),
            ','.join(str(ratio).split('.')),','.join(str(h_f).split('.'))),
        description='', type=ANALYSIS, atTime=None, waitMinutes=0, waitHours=0,
        queue=None, memory=90, memoryUnits=PERCENTAGE,
        getMemoryFromAnalysis=True, explicitPrecision=SINGLE,
        nodalOutputPrecision=FULL, echoPrint=OFF, modelPrint=OFF,
        contactPrint=OFF, historyPrint=OFF, userSubroutine='', scratch='',
        resultsFormat=ODB, multiprocessingMode=DEFAULT, numCpus=2,
        numDomains=2, numGPUs=0)

    # ------------------------------------------------------------------------
    # Save model

    # mdb.saveAs(pathName='A-{}R-{}H-{}'.format(','.join(str(h_u).split('.')),
    #         ','.join(str(ratio).split('.')),','.join(str(h_f).split('.'))))

if __name__ == '__main__':
    # createUndulation(h_u, h_f, ratio)
    createUndulation(0.2,0.2,0.2)

