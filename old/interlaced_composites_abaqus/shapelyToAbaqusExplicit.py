from sys import path
path.append(r'C:\Python27\Lib\site-packages')
path.append(r"C:\Users\rutge\Documents\GitHub\interlaced_model_creation")
from abaqus import *
from abaqusConstants import *
import regionToolset
import sigc as sigc
import tapePlacement as tp

session.viewports['Viewport: 1'].setValues(displayedObject=None)

# IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
# Define model parameters

cpt = 0.175  # cured ply thickness
laminateAngles = (0, 90, 45)  # define angles of tapes in laminate
tw = (25, 25, 25)  # tape widths
uw = 1.0

partTypes = ['Tape', 'Resin', 'Undulation']
# create grid using shapely
partGrid = sigc.createGrids(tapeAngles=laminateAngles, tapeWidths=tw,
                            undulationWidth=uw)
# identify parts in grid
tapePaths = tp.laminateCreation(
    grid=partGrid, tapeAngles=laminateAngles, tapeWidths=tw,
    tapeSpacing=1, undulationWidth=uw)

# IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
# Define model name

modelName = 'Interlaced Laminate'
mdb.Model(name=modelName, modelType=STANDARD_EXPLICIT)
laminateModel = mdb.models[modelName]
laminateAssembly = laminateModel.rootAssembly

# IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
# Materials

tapeMaterial = laminateModel.Material(name='Tape')
tapeMaterial.Density(table=((1.59e-06, ), ))
tapeMaterial.Depvar(deleteVar=20, n=24)
tapeMaterial.UserMaterial(
    mechanicalConstants=(146.8, 11.6, 11.6, 0.289, 0.289, 0.298, 6.47, 6.47,
                         4.38, 2.61, 1.759, 0.055, 0.285, 0.105, 53.0, 0.1,
                         0.1, 0.00075, 0.0025, 0.0035))
tapeSection = laminateModel.HomogeneousSolidSection(
        name='Tape Section', material='Tape')

undMaterial = laminateModel.Material(name='Undulation')
undMaterial.Depvar(deleteVar=20, n=24)
undMaterial.UserMaterial(
    mechanicalConstants=(72.36428, 5.8058, 5.8058, 0.2692, 0.2692, 0.3515,
                         3.2288, 14.3609, 8.7849, 2.61, 1.759, 0.055, 0.285,
                         0.105, 53.0, 0.1, 0.1, 0.00075, 0.0025, 0.0035))
undMaterial.Density(table=((1.59e-06, ), ))
undulationSection = laminateModel.HomogeneousSolidSection(
        name='Undulation Section', material='Undulation')

# explicit
# create step
laminateModel.ExplicitDynamicsStep(name='Loading Step', previous='Initial',
        timePeriod=2.5, massScaling=((SEMI_AUTOMATIC, MODEL, THROUGHOUT_STEP,
        0.0, 1e-06, BELOW_MIN, 1, 0, 0.0, 0.0, 0, None), ))

# Modify field output request
laminateModel.fieldOutputRequests['F-Output-1'].setValues(
        numIntervals=100)


# IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
# Generate parts

# This function is used to create parts using the geometric info from Shapely
def definePart(obj):
    x0, y0 = obj.exterior.xy
    abqCoords = zip(x0, y0)

    profileSketch = laminateModel.ConstrainedSketch(
            name='Sketch {}'.format(obj.objectID), sheetSize=200.0)
    for ind in range(len(abqCoords)-1):
        profileSketch.Line(
                point1=(abqCoords[ind][0], abqCoords[ind][1]),
                point2=(abqCoords[ind+1][0], abqCoords[ind+1][1]))

    # extrude profile to create part
    objPart = laminateModel.Part(name='Part {}'.format(obj.objectID),
            dimensionality=THREE_D, type=DEFORMABLE_BODY)
    objPart.BaseSolidExtrude(sketch=profileSketch, depth=cpt)

    objCells = objPart.cells  # all cells within part (only 1 cell)
    objRegion = regionToolset.Region(cells=objCells[:])

    if obj.objectType == 'Resin' or obj.objectType == 'Tape':
        objPart.SectionAssignment(
            region=objRegion, sectionName='Tape Section',
            offset=0.0, offsetType=MIDDLE_SURFACE, offsetField='',
            thicknessAssignment=FROM_SECTION)
    else:
        objPart.SectionAssignment(
            region=objRegion, sectionName='Undulation Section',
            offset=0.0, offsetType=MIDDLE_SURFACE, offsetField='',
            thicknessAssignment=FROM_SECTION)

    objPart.MaterialOrientation(
        region=objRegion, orientationType=SYSTEM, axis=AXIS_3,
        localCsys=None, fieldName='', additionalRotationType=ROTATION_ANGLE,
        additionalRotationField='', angle=obj.angle[-1],
        stackDirection=STACK_3)

    # create the part instance
    instanceName = 'Instance {}'.format(obj.objectID)
    objInstance = laminateAssembly.Instance(
            name=instanceName, part=objPart, dependent=ON)
    laminateAssembly.translate(
            instanceList=(instanceName, ),
            vector=(0.0, 0.0, (obj.layer-1)*cpt))
    return instanceName


# ----------------------------------------------------------------------------
# In the next section the Shapely geometry objects are iterated over and
# definePart() is called to create the parts

for layerNumber in range(1, len(partGrid)+1):
    for objID, obj in partGrid[layerNumber].iteritems():
        if obj.geom_type == 'Polygon':
            definePart(obj)
        elif obj.geom_type == 'MultiPolygon':
            for subObj in obj:
                definePart(subObj)

# this is rather difficult to read but all it does is create a flat
# list of the objIDs of each polygon in a tapepath so that they can be
# merged together in the next part of the program
for i, p in enumerate(tapePaths):
    if p:
        r = []
        temp = tapePaths[i+1:]
        flatList = [item for sublist in temp for item in sublist]
        for objID in p:
            if objID in flatList:
                r.append(objID)
        tapePaths[i] = [x for x in p if x not in r]
    else:
        continue

# this section merges all the parts in a tapepath into one part
for (pathNumber, tapePath) in enumerate(tapePaths):
    if len(tapePath) > 1:
        layer, gridID = tapePath[0].split('-')
        pathAngle = partGrid[int(layer)][int(gridID)].angle[-1]
        instanceList = [laminateAssembly.instances['Instance {}'.format(name)]
                        for name in tapePath]
        laminateAssembly.InstanceFromBooleanMerge(
            name='Pass {}'.format(pathNumber), instances=instanceList,
            keepIntersections=ON, originalInstances=DELETE, mergeNodes=ALL,
            nodeMergingTolerance=1e-06, domain=BOTH)
        laminateModel.parts['Pass {}'.format(pathNumber)].seedPart(
            size=1.0, deviationFactor=0.1, minSizeFactor=0.1)
        laminateModel.parts['Pass {}'.format(pathNumber)].generateMesh()
        passRegion = regionToolset.Region(
            cells=laminateModel.parts['Pass {}'.format(pathNumber)].cells[:])
        laminateModel.parts['Pass {}'.format(pathNumber)].MaterialOrientation(
            region=passRegion, orientationType=SYSTEM, axis=AXIS_3,
            localCsys=None, fieldName='', additionalRotationType=ROTATION_ANGLE,
            additionalRotationField='', angle=pathAngle, stackDirection=STACK_3)

# IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
# Create cohesive zone interactions between faces of the assembly
# Calculate CZM parameters according to Turon 2007 (sort of...)

alpha = 50
E_33 = tapeMaterial.userMaterial.mechanicalConstants[2]
K = (alpha*E_33)/cpt

# define contact properties
laminateModel.ContactProperty('Tangential')
laminateModel.interactionProperties['Tangential'].TangentialBehavior(
    formulation=PENALTY, directionality=ISOTROPIC, slipRateDependency=OFF,
    pressureDependency=OFF, temperatureDependency=OFF, dependencies=0,
    table=((0.15, ), ), shearStressLimit=None, maximumElasticSlip=FRACTION,
    fraction=0.005, elasticSlipStiffness=None)
laminateModel.interactionProperties['Tangential'].NormalBehavior(
    pressureOverclosure=HARD, allowSeparation=ON,
    constraintEnforcementMethod=DEFAULT)
laminateModel.ContactProperty('Cohesive')
laminateModel.interactionProperties['Cohesive'].CohesiveBehavior(
    defaultPenalties=OFF, table=((K, K, K), ))
laminateModel.interactionProperties['Cohesive'].Damage(
    criterion=QUAD_TRACTION, initTable=((0.055, 0.09, 0.09), ),
    useEvolution=ON, evolutionType=ENERGY, useMixedMode=ON,
    mixedModeType=POWER_LAW, exponent=1.0,
    evolTable=((3.52e-4, 3.52e-4, 3.52e-4), ))

# determine contacts
laminateModel.contactDetection(defaultType=CONTACT,
    interactionProperty='Cohesive', nameEachSurfaceFound=OFF,
    createUnionOfMasterSurfaces=ON, createUnionOfSlaveSurfaces=ON,
    separationTolerance=cpt/20.0)

# convert individual contact interactions to a single interaction with
# each contact defined as an 'Individual Property Assignment'
laminateModel.ContactExp(name='GC', createStepName='Initial')
laminateModel.interactions['GC'].includedPairs.setValuesInStep(
    stepName='Initial', useAllstar=ON)
laminateModel.interactions['GC'].contactPropertyAssignments.appendInStep(
        stepName='Initial', assignments=((GLOBAL, SELF, 'Tangential'), ))
for interact in laminateModel.interactions.values():
    if interact.name == 'GC':
        continue
    else:
        masterName = interact.master[0]
        slaveName = interact.slave[0]
        masterSurf = laminateAssembly.surfaces[masterName]
        slaveSurf = laminateAssembly.surfaces[slaveName]
        laminateModel.interactions['GC'].contactPropertyAssignments.appendInStep(
            stepName='Initial',
            assignments=((masterSurf, slaveSurf, 'Cohesive'), ))
        del laminateModel.interactions['{}'.format(interact.name)]

# IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII
# Apply boundary conditions

# identify faces at top and bottom for coupling
tFaces = [inst.faces.getByBoundingBox(-60.0, 74.0, -1.0, 60.0, 76.0, 6.0)
          for inst in laminateAssembly.instances.values()
          if inst.faces.getByBoundingBox(-60.0, 74.0, -1.0, 60.0, 76.0, 6.0)]
tFacesSurface = laminateAssembly.Surface(
        side1Faces=tFaces, name='Top Faces')
bFaces = [inst.faces.getByBoundingBox(-60.0, -76.0, -1.0, 60.0, -74.0, 6.0)
          for inst in laminateAssembly.instances.values()
          if inst.faces.getByBoundingBox(-60.0, -76.0, -1.0, 60.0, -74.0, 6.0)]
bFacesSurface = laminateAssembly.Surface(
        side1Faces=bFaces, name='Bottom Faces')
symFaces = [inst.faces.getByBoundingBox(-60.0, -76.0, -1.0, 60.0, 76.0, 0.001)
            for inst in laminateAssembly.instances.values()
            if inst.faces.getByBoundingBox(-60.0, -76.0, -1.0, 60.0, 76.0, 0.001)]
symFacesSet = laminateAssembly.Set(
        faces=symFaces, name='Symmetry Faces')

# create reference points for coupling
rfPoint1Id = laminateAssembly.ReferencePoint(
        point=(0.0, 80.0, 2.0)).id
rfPoint1 = laminateAssembly.referencePoints[rfPoint1Id]
rfPoint1Region = laminateAssembly.Set(
        referencePoints=(rfPoint1,), name='Coupling Reference Point 1')
rfPoint2Id = laminateAssembly.ReferencePoint(
        point=(0.0, -80.0, 2.0)).id
rfPoint2 = laminateAssembly.referencePoints[rfPoint2Id]
rfPoint2Region = laminateAssembly.Set(
        referencePoints=(rfPoint2,), name='Coupling Reference Point 2')

laminateModel.Coupling(
            name='Top Coupling', controlPoint=rfPoint1Region,
            surface=tFacesSurface, influenceRadius=WHOLE_SURFACE,
            couplingType=KINEMATIC, localCsys=None, u1=ON, u2=ON, u3=ON,
            ur1=ON, ur2=ON, ur3=ON)

laminateModel.Coupling(
            name='Bottom Coupling', controlPoint=rfPoint2Region,
            surface=bFacesSurface, influenceRadius=WHOLE_SURFACE,
            couplingType=KINEMATIC, localCsys=None, u1=ON, u2=ON, u3=ON,
            ur1=ON, ur2=ON, ur3=ON)

# explicit
laminateModel.SmoothStepAmplitude(name='Smoothing Amplitude',
    timeSpan=STEP, data=((0.0, 0.0), (1e-05, 1.0)))
laminateModel.DisplacementBC(
    name='Bottom Surface BC', createStepName='Loading Step',
    region=rfPoint2Region, u1=0.0, u2=0.0, u3=0.0, ur1=0.0, ur2=0.0,
    ur3=0.0, amplitude=UNSET, fixed=OFF, distributionType=UNIFORM,
    fieldName='', localCsys=None)
laminateModel.VelocityBC(name='Top Surface BC', createStepName='Loading Step',
    region=rfPoint1Region, v1=UNSET, v2=0.02, v3=UNSET, vr1=UNSET, vr2=UNSET,
    vr3=UNSET, amplitude='Smoothing Amplitude', localCsys=None,
    distributionType=UNIFORM, fieldName='')
laminateModel.ZsymmBC(name='Symmetry', createStepName='Loading Step',
    region=symFacesSet, localCsys=None)

mdb.saveAs(pathName=modelName)
