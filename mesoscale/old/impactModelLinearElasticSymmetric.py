from abaqus import *
from abaqusConstants import *
import regionToolset
import mesh
import itertools
import numpy as np
import math

# initialize viewport and name model
session.viewports['Viewport: 1'].setValues(displayedObject=None)
modelName = 'DWT_QI'   

mdb.Model(name=modelName, modelType=STANDARD_EXPLICIT)
impactModel = mdb.models[modelName]
impactAssembly = impactModel.rootAssembly

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<< PARAMETERS >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

couponL = 150.0 # coupon length [mm]
couponW = 100.0 # coupon width [mm]

impactorMass = 5.0 # mass of whole impactor (not 1/4) [kg]

cpt = 0.213
nLayers = 24
t = nLayers*cpt
angles = [0,45,-45,90]
m = (nLayers/int(len(angles)))/2
angleList = angles*m + list(reversed(angles*m))

time = 5.0 # duration to simulate [ms]
outputIntervals = 50 # requested field output intervals

# <<<<<<<<<<<<<<<<<<<<<<<<<<<< CREATE PARTS >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# Create the impactor
# Sketch the impactor (to create part by revolution)
impactorSketch = impactModel.ConstrainedSketch(name='Impactor Sketch', 
        sheetSize=200.0)
impactorSketch.ConstructionLine(point1=(0.0, -100.0), point2=(0.0, 100.0))
impactorSketch.ArcByCenterEnds(center=(0.0, 0.0), point1=(8.0, 0.0),
        point2=(0.0, -8.0), direction=CLOCKWISE)
impactorSketch.Line(point1=(0.0, -8.0), point2=(0.0, 32.0))
impactorSketch.Line(point1=(0.0, 32.0), point2=(8.0, 32.0))
impactorSketch.Line(point1=(8.0, 32.0), point2=(8.0, 0.0))
impactorPart = impactModel.Part(name='Impactor', dimensionality=THREE_D, 
        type=DISCRETE_RIGID_SURFACE)
impactorPart.BaseSolidRevolve(sketch=impactorSketch, angle=90.0,
        flipRevolveDirection=ON)
impactorRFPointID = impactorPart.ReferencePoint(point=(0.0, 0.0, 0.0)).id
impactorRFPoint = impactorPart.referencePoints[impactorRFPointID]
impactorRFPointRegion = impactorPart.Set(referencePoints=(impactorRFPoint,),
        name='Impactor Reference Point')

# Create shell
impactorCells = impactorPart.cells
impactorPart.RemoveCells(cellList = impactorCells[0:1])

# Assign inertia/mass to impactor
impactorPart.engineeringFeatures.PointMassInertia(
        name='Impactor Mass-Inertia', region=impactorRFPointRegion,
        mass=impactorMass/4.0, alpha=0.0, composite=0.0)

# ------------------------------------------------------------------------
# Create the clamps
clampCoords = (50.0,43.75)
clampSketch = impactModel.ConstrainedSketch(name='Clamp Sketch', 
        sheetSize=200.0)
clampSketch.CircleByCenterPerimeter(center=clampCoords,
        point1=(clampCoords[0]+4.0, clampCoords[1]))
clampPart = impactModel.Part(name='Clamp', dimensionality=THREE_D, 
        type=DISCRETE_RIGID_SURFACE)
clampPart.BaseSolidExtrude(sketch=clampSketch, depth=4.0)
clampRFPointID = clampPart.ReferencePoint(
    point=(clampCoords[0], clampCoords[1], 0.0)).id
clampRFPoint = clampPart.referencePoints[clampRFPointID]
clampRFPointRegion = clampPart.Set(referencePoints=(clampRFPoint,),
        name='Clamp Reference Point')

# Create shell
clampCells = clampPart.cells
clampPart.RemoveCells(cellList = clampCells[0:1])

# Assign inertia/mass to clamps
clampPart.engineeringFeatures.PointMassInertia(
        name='Clamp Mass-Inertia', region=clampRFPointRegion, mass=1.0,
        alpha=0.0, composite=0.0)

# ------------------------------------------------------------------------
# Create the bottom plate
bPlateCoords = [(0.0,37.5),(0.0,75.0),(100.0,75.0),(100.0,0.0),(67.5,0.0),
                (67.5,37.5),(0.0,37.5)]
bPlateSketch = impactModel.ConstrainedSketch(name='Bottom Plate Sketch', 
        sheetSize=200.0)
for k in range(len(bPlateCoords)-1):
    coord1 = bPlateCoords[k]
    coord2 = bPlateCoords[k+1]
    bPlateSketch.Line(point1=coord1, point2=coord2)
bPlatePart = impactModel.Part(name='Bottom Plate', dimensionality=THREE_D, 
        type=DISCRETE_RIGID_SURFACE)
bPlatePart.BaseShell(sketch=bPlateSketch)
bPlateRFPointID = bPlatePart.ReferencePoint(point=(0.0, 0.0, 0.0)).id
bPlateRFPoint = bPlatePart.referencePoints[bPlateRFPointID]
bPlateRFPointRegion = bPlatePart.Set(referencePoints=(bPlateRFPoint,),
        name='Bottom Plate Reference Point')

# Do not need to assign intertia as plate is constrained in every DOF

# ------------------------------------------------------------------------
# Create the specimen

specimenSketch = impactModel.ConstrainedSketch(name='Specimen Sketch', 
    sheetSize=200.0)
specimenSketch.rectangle(point1=(0.0, 0.0), point2=(couponL/2.0, couponW/2.0))
specimenPart = impactModel.Part(name='Specimen', dimensionality=THREE_D, 
    type=DEFORMABLE_BODY)
specimenPart.BaseSolidExtrude(sketch=specimenSketch, depth=t)
specimenCells = specimenPart.cells

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<< PARTITIONING >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# Through thickness partition into plies
plyDatumPlaneIds = []
for k in range(1,nLayers):
    plyDatumPlaneIds.append(specimenPart.DatumPlaneByPrincipalPlane(
        principalPlane=XYPLANE, offset=cpt*k).id)
    specimenCells = specimenPart.cells
    specimenPart.PartitionCellByDatumPlane(
            datumPlane=specimenPart.datums[plyDatumPlaneIds[k-1]],
            cells=specimenCells)

# ------------------------------------------------------------------------
# Partition for meshing

# Circular partition

# define the partition sketch face
partition_sketch_face = specimenPart.faces.findAt((25.0, 37.5, t), )
# choose an edge to define the orientation of the sketch partition
# "choose an edge vertically and to the right"
partition_edge = specimenPart.edges.findAt((75.0, 37.5, t), )
partition_definition = specimenPart.MakeSketchTransform(
        sketchPlane=partition_sketch_face, sketchUpEdge=partition_edge,
        sketchPlaneSide=SIDE1, origin=(0.0, 0.0, t))
partition_sketch = impactModel.ConstrainedSketch(
        name='Partition Sketch', sheetSize=4.54, gridSpacing=0.11,
        transform=partition_definition)
partition_sketch.rectangle(point1=(0.0, 0.0), point2=(20.0, 20.0))
specimenPart.projectReferencesOntoSketch(
            sketch=partition_sketch, filter=COPLANAR_EDGES)
specimenPart.PartitionCellBySketch(
            sketchPlane=partition_sketch_face, sketchUpEdge=partition_edge, 
            cells=specimenPart.cells, sketch=partition_sketch)
partition_path =  (specimenPart.edges.findAt(((20.0, 10.0, t), ),
            ((10.0, 20.0, t), )))
specimenPart.PartitionCellByExtrudeEdge(
    line=specimenPart.edges.findAt(coordinates=(75.0, 50.0, t-cpt/2.0)),
    cells=specimenPart.cells, edges=partition_path, sense=REVERSE)

# define the partition sketch face
partition_sketch_face = specimenPart.faces.findAt((25.0, 37.5, t), )
# choose an edge to define the orientation of the sketch partition
# "choose an edge vertically and to the right"
partition_edge = specimenPart.edges.findAt((75.0, 37.5, t), )
partition_definition = specimenPart.MakeSketchTransform(
        sketchPlane=partition_sketch_face, sketchUpEdge=partition_edge,
        sketchPlaneSide=SIDE1, origin=(0.0, 0.0, t))
partition_sketch2 = impactModel.ConstrainedSketch(
        name='Partition Sketch', sheetSize=4.54, gridSpacing=0.11,
        transform=partition_definition)
partition_sketch2.rectangle(point1=(0.0, 0.0), point2=(40.0, 40.0))
specimenPart.projectReferencesOntoSketch(
            sketch=partition_sketch2, filter=COPLANAR_EDGES)
specimenPart.PartitionCellBySketch(
            sketchPlane=partition_sketch_face, sketchUpEdge=partition_edge, 
            cells=specimenPart.cells, sketch=partition_sketch2)
partition_path2 =  (specimenPart.edges.findAt(((40.0, 20.0, t), ),
            ((20.0, 40.0, t), )))
specimenPart.PartitionCellByExtrudeEdge(
    line=specimenPart.edges.findAt(coordinates=(75.0, 50.0, t-cpt/2.0)),
    cells=specimenPart.cells, edges=partition_path2, sense=REVERSE)

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< MATERIALS >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# Create and assign material for the specimen
impactModel.Material(name='Material-1')
impactModel.materials['Material-1'].Density(table=((1.53e-06, ), ))
impactModel.materials['Material-1'].Elastic(
        type=ENGINEERING_CONSTANTS, table=((178.53, 8.10, 8.10, 0.252, 0.252, 0.392, 
        4.978, 4.978, 3.255), )) # estimated using Helius
impactModel.HomogeneousSolidSection(name='Section-1', 
    material='Material-1', thickness=None)
region = specimenPart.Set(cells=specimenCells, name='Specimen Cells')
specimenPart.SectionAssignment(region=region, sectionName='Section-1',
        offset=0.0, offsetType=MIDDLE_SURFACE, offsetField='', 
        thicknessAssignment=FROM_SECTION)

# Assign material orientations
for j in range(1,nLayers+1):
    layerAngle = angleList[j-1]
    layerCell = specimenPart.cells.getByBoundingBox(
            -76.0, -51.0, cpt*(j-1)-0.001, 76.0, 51.0, cpt*j+0.001)
    layerCellSet = specimenPart.Set(cells=layerCell, name='Layer {}'.format(j))
    specimenPart.MaterialOrientation(region=layerCellSet, 
        orientationType=SYSTEM, axis=AXIS_3, localCsys=None, 
        fieldName='', additionalRotationType=ROTATION_ANGLE, 
        additionalRotationField='', angle=layerAngle, stackDirection=STACK_3)

# <<<<<<<<<<<<<<<<<<<<<<<<<<< CREATE ASSEMBLY >>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# Create the part instances
impactAssembly = impactModel.rootAssembly
impactorInstance = impactAssembly.Instance(
        name='Impactor Instance', part=impactorPart, dependent=ON)
bPlateInstance = impactAssembly.Instance(
        name='Bottom Plate Instance', part=bPlatePart, dependent=ON)
specimenInstance = impactAssembly.Instance(
        name='Specimen Instance', part=specimenPart, dependent=ON)
clampInstance1 = impactAssembly.Instance(
        name='Clamp Instance 1', part=clampPart, dependent=ON)

impactAssembly.rotate(instanceList=('Bottom Plate Instance',
    'Specimen Instance', 'Clamp Instance 1', 'Clamp Instance 2',
    'Clamp Instance 3', 'Clamp Instance 4'), axisPoint=(0.0, 0.0, 0.0),
    axisDirection=(1.0, 0.0, 0.0), angle=-90.0)
impactAssembly.translate(instanceList=('Impactor Instance', ),
        vector=(0.0, 8.75+t, 0.0))
impactAssembly.translate(instanceList=('Clamp Instance 1', ), 
        vector=(0.0, t, 0.0))

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<<< SYMMETRY >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

impactorZSymmFace = impactorInstance.faces.findAt(((4.0, 35.0, 0.0), ))
specimenZSymmFaces = specimenInstance.faces.getByBoundingBox(
            -1.0, -1.0, -1.0, 76.0, t+1.0, 0.001)
bPlateZSymmEdge = bPlateInstance.edges.findAt(((75.0, 0.0, 0.0), ))
zSymmRegion = impactAssembly.Set(edges=bPlateZSymmEdge,
    faces=impactorZSymmFace+specimenZSymmFaces, name='Z-Symmetry')
impactModel.ZsymmBC(name='Z-Symmetry', createStepName='Initial', 
    region=zSymmRegion, localCsys=None)

impactorXSymmFace = impactorInstance.faces.findAt(((0.0, 35.0, -4.0), ))
specimenXSymmFaces = specimenInstance.faces.getByBoundingBox(
            -1.0, -1.0, -67.5, 0.001, t+1.0, 1.0)
bPlateXSymmEdge = bPlateInstance.edges.findAt(((0.0, 0.0, -67.5), ))
xSymmRegion = impactAssembly.Set(edges=bPlateXSymmEdge,
    faces=impactorXSymmFace+specimenXSymmFaces, name='X-Symmetry')
impactModel.XsymmBC(name='X-Symmetry', createStepName='Initial', 
    region=xSymmRegion, localCsys=None)

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<< CREATE STEP >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# create implicit step 
impactModel.ExplicitDynamicsStep(name='Loading Step',
        previous='Initial', timePeriod=time,
        massScaling=((SEMI_AUTOMATIC, MODEL, THROUGHOUT_STEP, 0.0,
            1e-06, BELOW_MIN, 1, 0, 0.0, 0.0, 0, None), ))

impactModel.fieldOutputRequests['F-Output-1'].setValues(
    variables=('S', 'LE', 'U', 'V', 'A', 'RF', 'CSTRESS',
        'CSDMG'), numIntervals=outputIntervals)

# <<<<<<<<<<<<<<<<<<<<<<<<<< APPLY CONSTRAINTS >>>>>>>>>>>>>>>>>>>>>>>>>>>>

# Encastre bottom plate
bPlateBCRegion = bPlateInstance.sets['Bottom Plate Reference Point']
impactModel.EncastreBC(name='Bottom Plate Encastre BC', 
        createStepName='Loading Step', region=bPlateBCRegion, localCsys=None)

# Fix clamps in all but y-direction
clampBCRFPoint1 = clampInstance1.referencePoints.values()[0]
clampBCRegion = impactAssembly.Set(
    referencePoints=(clampBCRFPoint1,), name='Clamp BC Region')
impactModel.DisplacementBC(name='Clamp Displacement BC', 
    createStepName='Loading Step', region=clampBCRegion, u1=0.0, u2=UNSET,
    u3=0.0, ur1=0.0, ur2=0.0, ur3=0.0, amplitude=UNSET, fixed=OFF, 
    distributionType=UNIFORM, fieldName='', localCsys=None)

# Apply load to clamps
impactModel.SmoothStepAmplitude(name='Smoothing Amplitude', timeSpan=STEP, 
    data=((0.0, 0.0), (1e-05, 1.0)))
impactModel.ConcentratedForce(name='Clamp Load BC', 
    createStepName='Loading Step', region=clampBCRegion, cf2=-1.0, 
    amplitude='Smoothing Amplitude', distributionType=UNIFORM, field='',
    localCsys=None)

# Fix impactor in all but y direction
impactorBCRegion = impactorInstance.sets['Impactor Reference Point']
impactModel.DisplacementBC(name='Impactor Displacement BC', 
    createStepName='Loading Step', region=impactorBCRegion, u1=0.0, u2=UNSET,
    u3=0.0, ur1=0.0, ur2=0.0, ur3=0.0, amplitude=UNSET, fixed=OFF, 
    distributionType=UNIFORM, fieldName='', localCsys=None)

# Apply velocity BC to impactor
impactModel.Velocity(name='ImpactorVelocity', region=impactorBCRegion, field='', 
    distributionType=MAGNITUDE, velocity1=0.0, velocity2=-3.4641, 
    velocity3=0.0, omega=0.0)

# <<<<<<<<<<<<<<<<<<<<<<<<< CREATE INTERACTIONS >>>>>>>>>>>>>>>>>>>>>>>>>>>

# Define contact properties
impactModel.ContactProperty('Contact Properties')
impactModel.interactionProperties['Contact Properties'].NormalBehavior(
    pressureOverclosure=HARD, allowSeparation=ON, 
    constraintEnforcementMethod=DEFAULT)
impactModel.interactionProperties['Contact Properties'].TangentialBehavior(
    formulation=PENALTY, directionality=ISOTROPIC, slipRateDependency=OFF, 
    pressureDependency=OFF, temperatureDependency=OFF, dependencies=0, 
    table=((0.15, ), ), shearStressLimit=None, maximumElasticSlip=FRACTION, 
    fraction=0.005, elasticSlipStiffness=None)
# Define general contact
impactModel.ContactExp(name='GC', 
    createStepName='Loading Step')
impactModel.interactions['GC'].includedPairs.setValuesInStep(
    stepName='Loading Step', useAllstar=ON)
impactModel.interactions['GC'].contactPropertyAssignments.appendInStep(
    stepName='Loading Step', assignments=((GLOBAL, SELF, 'Contact Properties'), ))

# <<<<<<<<<<<<<<<<<<<<<<<<<<<< CREATE MESH >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

# Mesh impactor 
elemType1 = mesh.ElemType(elemCode=R3D4, elemLibrary=STANDARD)
elemType2 = mesh.ElemType(elemCode=R3D3, elemLibrary=STANDARD)
impactorFaces = impactorPart.faces
pickedRegions =(impactorFaces, )
impactorPart.setElementType(regions=pickedRegions,
        elemTypes=(elemType1, elemType2))
impactorPart.seedPart(size=1.0, deviationFactor=0.1, minSizeFactor=0.1)
impactorPart.generateMesh()

# Mesh clamp 
clampFaces = clampPart.faces
pickedRegions =(clampFaces, )
clampPart.setElementType(regions=pickedRegions,
        elemTypes=(elemType1, elemType2))
clampPart.seedPart(size=1.0, deviationFactor=0.1, minSizeFactor=0.1)
clampPart.generateMesh()

# Mesh bottom plate
bPlateFaces = bPlatePart.faces
pickedRegions =(bPlateFaces, )
bPlatePart.setElementType(regions=pickedRegions,
        elemTypes=(elemType1, elemType2))
bPlatePart.seedPart(size=2.5, deviationFactor=0.1, minSizeFactor=0.1)
bPlatePart.generateMesh()

# ------------------------------------------------------------------------
# Mesh specimen 

# Mesh sizes
nOutside = 5.0
nMiddle = 1.0
nTarget = 0.3

# First assign global mesh seeds to seed the outside edges
specimenPart.seedPart(size=nOutside, deviationFactor=0.1, minSizeFactor=0.1)

for p in range(nLayers+1):
    middleEdges = specimenPart.edges.findAt(((30.0, 0.0, cpt*p), ),
            ((0.0, 30.0, cpt*p), ), ((40.0, 20.0, cpt*p), ), ((20.0, 40.0, cpt*p), ))
    specimenPart.seedEdgeBySize(edges=middleEdges, size=nMiddle, constraint=FINER)
    targetEdges = specimenPart.edges.findAt(((10.0, 0.0, cpt*p), ),
            ((0.0, 10.0, cpt*p), ), ((10.0, 20.0, cpt*p), ), ((20.0, 10.0, cpt*p), ))
    specimenPart.seedEdgeBySize(edges=targetEdges, size=nTarget, constraint=FINER)

specimenPart.generateMesh()

# <<<<<<<<<<<<<<<<<<<<<<<<<<<<< CREATE JOB >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Create job

mdb.Job(name=modelName, model=modelName, description='', type=ANALYSIS,
        atTime=None, waitMinutes=0, waitHours=0, queue=None, memory=90, 
        memoryUnits=PERCENTAGE, explicitPrecision=DOUBLE_PLUS_PACK, 
        nodalOutputPrecision=FULL, echoPrint=OFF, modelPrint=OFF, 
        contactPrint=OFF, historyPrint=OFF, userSubroutine='', scratch='', 
        resultsFormat=ODB, parallelizationMethodExplicit=DOMAIN, numDomains=4, 
        activateLoadBalancing=False, multiprocessingMode=DEFAULT, numCpus=4)

# Save model

# mdb.saveAs(pathName=modelName)
