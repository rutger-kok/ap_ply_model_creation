from abaqus import *
from abaqusConstants import *
from caeModules import *
from visualization import *
from math import pi
import mesh
import regionToolset
import os
from sys import path
import numpy as np
path.append('C:\\Python27\\Lib\\site-packages')
path.append('C:\\Workspace\\ELSA')
from shapely.geometry import LineString
from shapely.geometry import Point

class TensileModel():
    def __init__(self,t, nPlies, Lo, Lg, D, wo, wg):
        '''
        Set specimen parameters
        Lo = overall specimen length (mm)
        Lg = length of gauge section (mm)
        D = distance between end tabs (mm)
        wo = overall specimen width (mm)
        wg = width of gauge section (mm)
        r1 = radius of shoulder 1 (mm)
        '''
        self.t = t
        self.nPlies = nPlies
        self.Lo = Lo
        self.Lg = Lg
        self.D = D
        self.wo = wo
        self.wg = wg
        self.meshSize = 1.0
        self.r = ((self.D**2 - 2 * self.D * self.Lg + self.wo**2 -
                  2 * self.wo * self.wg + self.Lg**2 + self.wg**2) /
                  (4 * (self.wo - self.wg)))

        # Set name
        self.modelName = 'Lg-{}_wo-{}_r-{}'.format(int(Lg), int(wo), int(self.r))
        mdb.Model(name=self.modelName, modelType=STANDARD_EXPLICIT)
        self.model = mdb.models[self.modelName]
        self.assembly = self.model.rootAssembly
        # set the work directory
        wd = 'C:\\Workspace\\ELSA\\{}'.format(self.modelName)
        if not os.path.exists(wd):
            os.makedirs(wd)
        os.chdir(wd)


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

        mat = self.model.Material(name='VTC-401')
        mat.Elastic(type=LAMINA, table=((E11, E22, nu12, G12, G13, G23), ))


    def createSpecimenPart(self):
        partSketch = self.model.ConstrainedSketch(name='__profile__', 
            sheetSize=200.0)
        g = partSketch.geometry
        # leftmost edge
        partSketch.Line(point1=(0.0, self.wo / 2.0), point2=(0.0, self.wo))
        # top left
        partSketch.Line(point1=(0.0, self.wo),
                        point2=((self.Lo - self.D) / 2.0, self.wo))
        # bottom line
        partSketch.Line(point1=(0.0, self.wo / 2.0),
                        point2=(self.Lo, self.wo / 2.0))
        # rightmost edge
        partSketch.Line(point1=(self.Lo, self.wo / 2.0), point2=(self.Lo, self.wo))
        # top right 
        partSketch.Line(point1=(self.Lo, self.wo),
                        point2=((self.Lo + self.D) / 2.0, self.wo))
        # gauge length line
        gl = partSketch.Line(point1=((self.Lo - self.Lg) / 2.0,
                                (self.wo / 2.0) + (self.wg / 2.0)),
                             point2=((self.Lo + self.Lg) / 2.0,
                                (self.wo / 2.0) + (self.wg / 2.0)))

        arc1 = partSketch.ArcByStartEndTangent(
            point1=((self.Lo - self.Lg) / 2.0,
                    (self.wo / 2.0) + (self.wg / 2.0)),
            point2=((self.Lo - self.D) / 2.0, self.wo),
            entity=gl)

        mirrorLine1 = partSketch.ConstructionLine(point1=(self.Lo / 2.0, 0.0),
                                    point2=(self.Lo /2.0, 100.0))
        partSketch.copyMirror(
            mirrorLine=g.findAt((self.Lo / 2.0, self.wo  / 2.0 + self.wg / 4.0)),
            objectList=(g.findAt(arc1.pointOn), ))                        

        part = self.model.Part(name='Specimen', dimensionality=THREE_D,
            type=DEFORMABLE_BODY)
        part.BaseShell(sketch=partSketch)

        # partition part into clamped and unclamped regions
        dplaneID1 = part.DatumPlaneByPrincipalPlane(
            principalPlane=YZPLANE, offset=(self.Lo - self.D) / 2.0).id
        part.PartitionFaceByDatumPlane(
            datumPlane=part.datums[dplaneID1], faces=part.faces)
        dplaneID2 = part.DatumPlaneByPrincipalPlane(
            principalPlane=YZPLANE, offset=(self.Lo + self.D) / 2.0).id
        part.PartitionFaceByDatumPlane(
            datumPlane=part.datums[dplaneID2], faces=part.faces)

        # partition gauge section
        dplaneID3 = part.DatumPlaneByPrincipalPlane(
            principalPlane=YZPLANE, offset=(self.Lo - self.Lg) / 2.0).id
        part.PartitionFaceByDatumPlane(
            datumPlane=part.datums[dplaneID3], faces=part.faces)
        dplaneID4 = part.DatumPlaneByPrincipalPlane(
            principalPlane=YZPLANE, offset=(self.Lo + self.Lg) / 2.0).id
        part.PartitionFaceByDatumPlane(
            datumPlane=part.datums[dplaneID4], faces=part.faces)  

        faces = part.faces  # all faces within part
        faceRegion = regionToolset.Region(faces=faces[:])
        compositeLayup = part.CompositeLayup(
            name='XP20', description='', elementType=SHELL, 
            offsetType=MIDDLE_SURFACE, symmetric=True, 
            thicknessAssignment=FROM_SECTION)
        compositeLayup.Section(preIntegrate=OFF, integrationRule=SIMPSON, 
            thicknessType=UNIFORM, poissonDefinition=DEFAULT, temperature=GRADIENT, 
            useDensity=OFF)
        compositeLayup.ReferenceOrientation(orientationType=GLOBAL, localCsys=None, 
            fieldName='', additionalRotationType=ROTATION_NONE, angle=0.0, 
            axis=AXIS_3)
        for k in range(self.nPlies / 2):
            if k % 2 != 0:
                orient = 90.0
            else:
                orient = 0.0
            compositeLayup.CompositePly(suppressed=False, plyName='Ply-{}'.format(k+1), 
                region=faceRegion, material='VTC-401', thicknessType=SPECIFY_THICKNESS, 
                thickness=self.t, orientationType=SPECIFY_ORIENT, orientationValue=orient, 
                additionalRotationType=ROTATION_NONE, additionalRotationField='', 
                axis=AXIS_3, angle=0.0, numIntPoints=3)

        # mesh the part
        elemType1 = mesh.ElemType(
            elemCode=S4R, elemLibrary=STANDARD, secondOrderAccuracy=OFF,
            hourglassControl=ENHANCED)
        elemType2 = mesh.ElemType(elemCode=S3, elemLibrary=STANDARD)
        part.setElementType(regions=faceRegion, elemTypes=(elemType1, elemType2))
        part.seedPart(size=self.meshSize)
        part.generateMesh()

        # create the part instance
        self.assembly.Instance(name='Specimen-Instance', part=part, dependent=ON)
        self.specInst = self.assembly.instances['Specimen-Instance']

        # create gauge section set
        gaugeFace = self.specInst.faces.findAt(
            ((self.Lo / 2.0, self.wo / 2.0 + self.wg / 4.0, 0.0), ))
        gaugeFaceSet = self.assembly.Set(faces=gaugeFace, name='Gauge Section')

        betweenClamps = self.specInst.faces.getByBoundingBox(
            (self.Lo - self.D) / 2.0 - 0.01, -0.01, -0.01,
            (self.Lo + self.D) / 2.0 + 0.01, self.wo + 0.01, 0.01)
        betweenClampsSet = self.assembly.Set(faces=betweenClamps, name='Between Clamps')

    def createImplicitStep(self):
        # mass scaling applied if increment < 1e-6
        self.model.StaticStep(name='Loading Step', previous='Initial')
        self.model.fieldOutputRequests['F-Output-1'].setValues(
            variables=('S', 'E', 'LE', 'U', 'RF', 'CF'))


    def applyContraints(self):

        # identify edges at right for coupling
        rightFace = self.specInst.faces.findAt(((self.Lo / 2.0 + 3* self.D / 4.0,
                                                    self.wo / 2.0, 0.0), ))
        rightFaceSurface = self.assembly.Surface(side1Faces=rightFace,
                                                 name='Right Face')

        leftFace = self.specInst.faces.findAt(((self.Lo / 2.0 - 3* self.D / 4.0,
                                                   self.wo / 2.0, 0.0), ))
        leftFaceSurface = self.assembly.Surface(side1Faces=leftFace,
                                                 name='Left Face')

        symEdges = self.specInst.edges.getByBoundingBox(
            -0.01, -0.01, -0.01, self.Lo + 0.01, (self.wo / 2.0) + 0.01, 0.01)
        symEdgesSet = self.assembly.Set(edges=symEdges, name='Symmetry Edges')

        # create reference points for coupling
        rfPoint1Id = self.assembly.ReferencePoint(
            point=(-1.0, self.wo / 2.0, 0.0)).id
        rfPoint1 = self.assembly.referencePoints[rfPoint1Id]
        rfPoint1Region = self.assembly.Set(referencePoints=(rfPoint1,),
                                           name='Coupling Reference Point 1')
        rfPoint2Id = self.assembly.ReferencePoint(
            point=(self.Lo + 1.0, self.wo / 2.0, 0.0)).id
        rfPoint2 = self.assembly.referencePoints[rfPoint2Id]
        rfPoint2Region = self.assembly.Set(referencePoints=(rfPoint2,),
                                           name='Coupling Reference Point 2')

        self.model.Coupling(
            name='Right Coupling', controlPoint=rfPoint2Region,
            surface=rightFaceSurface, influenceRadius=WHOLE_SURFACE,
            couplingType=KINEMATIC, localCsys=None, u1=ON, u2=ON, u3=ON,
            ur1=ON, ur2=ON, ur3=ON)

        self.model.Coupling(
            name='Left Coupling', controlPoint=rfPoint1Region,
            surface=leftFaceSurface, influenceRadius=WHOLE_SURFACE,
            couplingType=KINEMATIC, localCsys=None, u1=ON, u2=ON,
            u3=ON, ur1=ON, ur2=ON, ur3=ON)

        self.model.SmoothStepAmplitude(name='Smoothing Amplitude', timeSpan=STEP,
                                       data=((0.0, 0.0), (1e-05, 1.0)))
        self.model.DisplacementBC(
            name='Clamping BC', createStepName='Loading Step',
            region=rfPoint1Region, u1=0.0, u2=0.0, u3=0.0, ur1=0.0, ur2=0.0,
            ur3=0.0, amplitude=UNSET, fixed=OFF, distributionType=UNIFORM,
            fieldName='', localCsys=None)
        self.model.DisplacementBC(
            name='Displacement BC', createStepName='Loading Step',
            region=rfPoint2Region, u1=1.0, u2=0.0, u3=0.0, ur1=0.0, ur2=0.0,
            ur3=0.0, amplitude=UNSET, fixed=OFF, distributionType=UNIFORM,
            fieldName='', localCsys=None)
        self.model.YsymmBC(name='Symmetry', createStepName='Loading Step',
                           region=symEdgesSet, localCsys=None)

   
    def createJob(self, cpus=4, run=False):
        mdb.Job(
            name=self.modelName, model=self.modelName, description='',
            atTime=None, waitMinutes=0, waitHours=0, queue=None, type=ANALYSIS,
            memoryUnits=PERCENTAGE, explicitPrecision=DOUBLE_PLUS_PACK,
            nodalOutputPrecision=FULL, echoPrint=OFF, modelPrint=OFF,
            contactPrint=OFF, historyPrint=OFF, userSubroutine='', scratch='',
            resultsFormat=ODB, parallelizationMethodExplicit=DOMAIN,
            numDomains=cpus, activateLoadBalancing=False,
            multiprocessingMode=DEFAULT, numCpus=cpus, memory=90)
        if run is True:
            mdb.jobs[self.modelName].submit(consistencyChecking=OFF)
            mdb.jobs[self.modelName].waitForCompletion()

    def saveModel(self):
        mdb.saveAs(pathName=self.modelName)

    
    def postProcess(self):

        # Open the Output Database for the current Job
        odb = openOdb(path='{}.odb'.format(self.modelName))

        frame = odb.steps['Loading Step'].frames[-1]
        sFieldOutput = frame.fieldOutputs['S']

        # determine nominal stress
        gaugeRegion = odb.rootAssembly.elementSets['Gauge Section']
        s_gaugeRegion = sFieldOutput.getSubset(region=gaugeRegion)
        s_gaugeRegion_mises = [x.mises for x in s_gaugeRegion.values]
        s_avg = np.mean(s_gaugeRegion_mises)
        # determine max stress between clamps
        betweenClamps = odb.rootAssembly.elementSets['Between Clamps']
        s_betweenClamps = sFieldOutput.getSubset(region=betweenClamps)
        s_betweenClamps_mises = [x.mises for x in s_betweenClamps.values]
        s_max = np.max(s_betweenClamps_mises)
        # calculate stress intensity factor
        sif = s_max / s_avg
        return s_max, sif


def main(plyThickness, nLayers, Lo, Lg, D, wo, wg):
    # Create model
    mdl = TensileModel(plyThickness, nLayers, Lo, Lg, D, wo, wg)
    mdl.createSpecimenPart()
    mdl.createMaterials()
    mdl.createImplicitStep()
    mdl.applyContraints()
    mdl.createJob(run=True)
    mdl.saveModel()
    smax, sif = mdl.postProcess()
    return smax, sif, int(mdl.r)

if __name__ == '__main__':
    smax, sif, r = main(0.2, 20, 300.0, 50.0, 100.0, 50.0, 40.0)
    print smax, sif, r


    
