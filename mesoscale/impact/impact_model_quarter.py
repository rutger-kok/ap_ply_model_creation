'''
Module used to define the experimental setup for a drop weight tower test,
i.e. the definition of the impactor, clamps, bottom panel etc.

(c) Rutger Kok, 27/11/2020
'''
from abaqus import *
from abaqusConstants import *
import mesh
from sys import path
githubPath = 'C:\\GitHub'
path.append('C:\\Python27\\Lib\\site-packages')
path.append(githubPath + '\\interlaced_model_creation\\mesoscale')
from interlaced_3d import Interlaced3D
from interlaced_2d import Interlaced2D


class ImpactModel(Interlaced3D, Interlaced2D):
    def __init__(self, model_name, dir_name=None):
        Interlaced2D.__init__(self, model_name, dir_name)

    def set_test_parameters(self, time, output_intervals, energy,
                            coarse_mesh=1.0, medium_mesh=0.5, fine_mesh=0.4):
        '''Set drop weight tower test parameters'''
        self.time = time
        self.output_intervals = output_intervals
        self.energy = energy
        self.impactor_mass = 0.00585  # tons
        # calculate initial impactor velocity
        self.init_velocity = -3280.0  # -((2.0 * self.energy) / self.impactor_mass)**0.5
        # mesh densities
        self.coarse_mesh = coarse_mesh
        self.medium_mesh = medium_mesh
        self.fine_mesh = fine_mesh

    def create_impactor_part(self):
        ''' Creates the hemispherical impactor for DWT simulations'''
        # Sketch the impactor (to create part by revolution)
        sketch = self.model.ConstrainedSketch(
            name='Impactor Sketch', sheetSize=200.0)
        sketch.ConstructionLine(point1=(0.0, -100.0), point2=(0.0, 100.0))
        sketch.ArcByCenterEnds(center=(0.0, 0.0), point1=(8.0, 0.0),
                               point2=(0.0, -8.0), direction=CLOCKWISE)
        sketch.Line(point1=(0.0, -8.0), point2=(0.0, 32.0))
        sketch.Line(point1=(0.0, 32.0), point2=(8.0, 32.0))
        sketch.Line(point1=(8.0, 32.0), point2=(8.0, 0.0))
        self.impactor_part = self.model.Part(
            name='Impactor', dimensionality=THREE_D,
            type=DISCRETE_RIGID_SURFACE)
        self.impactor_part.BaseSolidRevolve(
            sketch=sketch, angle=90.0, flipRevolveDirection=ON)
        rf_point_id = self.impactor_part.ReferencePoint(point=(0, 0, 0)).id
        rf_point = self.impactor_part.referencePoints[rf_point_id]
        rf_point_region = self.impactor_part.Set(
            referencePoints=(rf_point,), name='Impactor Reference Point')
        # Create shell
        cells = self.impactor_part.cells
        self.impactor_part.RemoveCells(cellList=cells[0:1])
        # Assign inertia/mass to impactor
        self.impactor_part.engineeringFeatures.PointMassInertia(
            name='Impactor Mass-Inertia', region=rf_point_region,
            mass=self.impactor_mass / 4.0, alpha=0.0, composite=0.0)
        # Meshing
        elem_type_1 = mesh.ElemType(elemCode=R3D4, elemLibrary=STANDARD)
        elem_type_2 = mesh.ElemType(elemCode=R3D3, elemLibrary=STANDARD)
        faces = self.impactor_part.faces
        picked_regions = (faces, )
        self.impactor_part.setElementType(
            regions=picked_regions, elemTypes=(elem_type_1, elem_type_2))
        self.impactor_part.seedPart(
            size=1.0, deviationFactor=0.1, minSizeFactor=0.1)
        self.impactor_part.generateMesh()

    def create_clamp_part(self):
        '''Create the restraining clamps for DWT simulations'''
        coords = (50.0, 44.0)
        sketch = self.model.ConstrainedSketch(name='Clamp Sketch',
                                              sheetSize=200.0)
        sketch.CircleByCenterPerimeter(center=coords,
                                       point1=(coords[0] + 4.0, coords[1]))
        self.clamp_part = self.model.Part(
            name='Clamp', dimensionality=THREE_D, type=DISCRETE_RIGID_SURFACE)
        self.clamp_part.BaseSolidExtrude(sketch=sketch, depth=4.0)
        rf_point_id = self.clamp_part.ReferencePoint(
            point=(coords[0], coords[1], 0.0)).id
        rf_point = self.clamp_part.referencePoints[rf_point_id]
        rf_point_region = self.clamp_part.Set(
            referencePoints=(rf_point,), name='Clamp Reference Point')
        # Create shell
        cells = self.clamp_part.cells
        self.clamp_part.RemoveCells(cellList=cells[0:1])
        # Assign inertia/mass to clamps
        self.clamp_part.engineeringFeatures.PointMassInertia(
            name='Clamp Mass-Inertia', region=rf_point_region, mass=0.001,
            alpha=0.0, composite=0.0)
        # Mesh clamp
        faces = self.clamp_part.faces
        picked_regions = (faces, )
        elem_type_1 = mesh.ElemType(elemCode=R3D4, elemLibrary=STANDARD)
        elem_type_2 = mesh.ElemType(elemCode=R3D3, elemLibrary=STANDARD)
        self.clamp_part.setElementType(
            regions=picked_regions, elemTypes=(elem_type_1, elem_type_2))
        self.clamp_part.seedPart(
            size=1.0, deviationFactor=0.1, minSizeFactor=0.1)
        self.clamp_part.generateMesh()

    def create_bottom_plate(self):
        '''Create the plate on which specimens rest for DWT simulations'''
        coords = [(0.0, 37.5), (0.0, 75.0), (100.0, 75.0), (100.0, 0.0),
                  (67.5, 0.0), (67.5, 37.5), (0.0, 37.5)]
        sketch = self.model.ConstrainedSketch(name='Bottom Plate Sketch',
                                              sheetSize=200.0)
        for k in range(len(coords) - 1):
            coord_1 = coords[k]
            coord_2 = coords[k + 1]
            sketch.Line(point1=coord_1, point2=coord_2)
        self.b_plate_part = self.model.Part(
            name='Bottom Plate', dimensionality=THREE_D,
            type=DISCRETE_RIGID_SURFACE)
        self.b_plate_part.BaseShell(sketch=sketch)
        rf_point_id = self.b_plate_part.ReferencePoint(point=(0, 0, 0)).id
        rf_point = self.b_plate_part.referencePoints[rf_point_id]
        self.b_plate_part.Set(referencePoints=(rf_point,),
                              name='Bottom Plate Reference Point')
        # Do not need to assign intertia as plate is constrained in every DOF
        # Mesh bottom plate
        faces = self.b_plate_part.faces
        picked_regions = (faces, )
        elem_type_1 = mesh.ElemType(elemCode=R3D4, elemLibrary=STANDARD)
        elem_type_2 = mesh.ElemType(elemCode=R3D3, elemLibrary=STANDARD)
        self.b_plate_part.setElementType(
            regions=picked_regions, elemTypes=(elem_type_1, elem_type_2))
        self.b_plate_part.seedPart(
            size=2.5, deviationFactor=0.1, minSizeFactor=0.1)
        self.b_plate_part.generateMesh()

    def create_assembly(self):
        # Create the part instances
        self.assembly.Instance(name='Impactor Instance',
                               part=self.impactor_part, dependent=ON)
        self.assembly.Instance(name='Bottom Plate Instance',
                               part=self.b_plate_part, dependent=ON)
        self.assembly.Instance(name='Clamp Instance', part=self.clamp_part,
                               dependent=ON)
        self.assembly.Instance(name='Specimen Instance',
                               part=self.specimen_part, dependent=ON)
        self.assembly.rotate(
            instanceList=('Bottom Plate Instance', 'Clamp Instance',
                          'Specimen Instance'),
            axisPoint=(0.0, 0.0, 0.0), axisDirection=(1.0, 0.0, 0.0),
            angle=-90.0)
        self.assembly.translate(instanceList=('Impactor Instance', ),
                                vector=(0.0, 8.75 + self.l_thickness, 0.0))
        self.assembly.translate(instanceList=('Clamp Instance', ),
                                vector=(0.0, self.l_thickness, 0.0))

    def constrain_b_plate(self):
        # Encastre bottom plate
        instance = self.assembly.instances['Bottom Plate Instance']
        bc_region = instance.sets['Bottom Plate Reference Point']
        self.model.EncastreBC(name='Bottom Plate Encastre BC',
                              createStepName='Loading Step', localCsys=None,
                              region=bc_region)

    def constrain_clamps(self):
        # Fix clamps in all but y-direction
        instance = self.assembly.instances['Clamp Instance']
        bc_region = instance.sets['Clamp Reference Point']
        self.model.DisplacementBC(
            name='Clamp Displacement BC', createStepName='Loading Step',
            region=bc_region, u1=0.0, u2=UNSET, u3=0.0, ur1=0.0, ur2=0.0,
            ur3=0.0, amplitude=UNSET, fixed=OFF, distributionType=UNIFORM,
            fieldName='', localCsys=None)
        # Apply load to clamps
        self.model.SmoothStepAmplitude(
            name='Smoothing Amplitude', timeSpan=STEP,
            data=((0.0, 0.0), (1e-05, 1.0)))
        self.model.ConcentratedForce(
            name='Clamp Load BC', createStepName='Loading Step',
            region=bc_region, cf2=-1100.0, amplitude='Smoothing Amplitude',
            distributionType=UNIFORM, field='', localCsys=None)

    def constrain_impactor(self):
        # Fix impactor in all but y direction
        impactor_instance = self.assembly.instances['Impactor Instance']
        impactor_bc_region = impactor_instance.sets['Impactor Reference Point']
        self.model.DisplacementBC(
            name='Impactor Displacement BC', createStepName='Loading Step',
            region=impactor_bc_region, u1=0.0, u2=UNSET, u3=0.0, ur1=0.0,
            ur2=0.0, ur3=0.0, amplitude=UNSET, fixed=OFF, fieldName='',
            distributionType=UNIFORM, localCsys=None)
        # Apply velocity BC to impactor
        self.model.Velocity(
            name='ImpactorVelocity', region=impactor_bc_region, field='',
            distributionType=MAGNITUDE, velocity1=0.0,
            velocity2=self.init_velocity, velocity3=0.0, omega=0.0)

    def create_test_specimen(self):
        sketch = self.model.ConstrainedSketch(
            name='Specimen Sketch', sheetSize=200.0)
        sketch.rectangle(point1=(0.0, 0.0), point2=(75.0, 50.0))
        part = self.model.Part(
            name='Specimen', dimensionality=THREE_D, type=DEFORMABLE_BODY)
        part.BaseSolidExtrude(sketch=sketch, depth=self.l_thickness)

        # assign properties
        self.model.HomogeneousSolidSection(
            name='Aluminium Section', material='Aluminium', thickness=None)
        cells = part.cells.findAt(((75.0, 33.333333, 1.186667), ))
        region = part.Set(cells=cells, name='Specimen Region')
        part.SectionAssignment(
            region=region, sectionName='Aluminium Section', offset=0.0,
            offsetType=MIDDLE_SURFACE, offsetField='',
            thicknessAssignment=FROM_SECTION)

        self.mesh_part_3d(part, self.fine_mesh)
        self.specimen_part = part

    def apply_symmetry(self):
        impactor_instance = self.assembly.instances['Impactor Instance']
        specimen_instance = self.assembly.instances['Specimen Instance']
        bplate_instance = self.assembly.instances['Bottom Plate Instance']
        impactor_z_sym = impactor_instance.faces.findAt(((4.0, 35.0, 0.0), ))
        specimen_z_sym = specimen_instance.faces.findAt(((4.0, 0.01, 0.0), ))
        bottom_plate_z_sym = bplate_instance.edges.findAt(((75.0, 0.0, 0.0), ))
        z_sym_region = self.assembly.Set(
            edges=bottom_plate_z_sym, faces=impactor_z_sym + specimen_z_sym,
            name='Z-Symmetry')
        self.model.ZsymmBC(
            name='Z-Symmetry', createStepName='Initial',
            region=z_sym_region, localCsys=None)

        impactor_x_sym = impactor_instance.faces.findAt(((0.0, 35.0, -4.0), ))
        specimen_x_sym = specimen_instance.faces.findAt(((0.0, 0.01, -4.0), ))
        bottom_plate_x_sym = bplate_instance.edges.findAt(
            ((0.0, 0.0, -67.5), ))
        x_sym_region = self.assembly.Set(
            edges=bottom_plate_x_sym, faces=impactor_x_sym + specimen_x_sym,
            name='X-Symmetry')
        self.model.XsymmBC(
            name='X-Symmetry', createStepName='Initial', region=x_sym_region,
            localCsys=None)

    def create_aluminium(self):
        aluminium = self.model.Material(name='Aluminium')
        aluminium.Density(table=((2.8e-09, ), ))
        aluminium.Elastic(table=((73858.0, 0.33), ))
        aluminium.Plastic(
            hardening=JOHNSON_COOK,
            table=((369.0, 684.0, 0.73, 1.0, 775.0, 293.0), ))
        aluminium.plastic.RateDependent(
            type=JOHNSON_COOK, table=((0.0083, 1.0), ))
        aluminium.JohnsonCookDamageInitiation(
            table=((0.112, 0.123, 0.5, 0.007, 0.0, 775.0, 273.0, 1.0), ))
        aluminium.johnsonCookDamageInitiation.DamageEvolution(
            type=DISPLACEMENT, table=((0.5, ), ))

    def create_interactions_explicit(self):
        interaction_props = self.model.ContactProperty('IntProp-1')
        interaction_props.TangentialBehavior(
            formulation=PENALTY, directionality=ISOTROPIC,
            slipRateDependency=OFF, pressureDependency=OFF,
            temperatureDependency=OFF, dependencies=0,
            table=((0.2, ), ), shearStressLimit=None, fraction=0.005,
            maximumElasticSlip=FRACTION, elasticSlipStiffness=None)
        interaction_props.NormalBehavior(
            pressureOverclosure=HARD, allowSeparation=ON,
            constraintEnforcementMethod=DEFAULT)
        general_contact = self.model.ContactExp(
            name='Int-1', createStepName='Initial')
        general_contact.includedPairs.setValuesInStep(
            stepName='Initial', useAllstar=ON)
        general_contact.contactPropertyAssignments.appendInStep(
            stepName='Initial',
            assignments=((GLOBAL, SELF, 'IntProp-1'), ))


if __name__ == '__main__':

    # Simulation parameters
    time = 0.025  # duration to simulate [s]
    output_intervals = 100  # requested field output intervals
    energy = 30000.0  # impact energy to simulate [Nmm]

    # Create model
    mdl = ImpactModel('DWT_Test_Quarter')
    mdl.set_test_parameters(time, output_intervals, energy)
    mdl.l_thickness = 1.78  # for test purposes
    mdl.create_aluminium()
    mdl.create_impactor_part()
    mdl.create_test_specimen()
    mdl.create_clamp_part()
    mdl.create_bottom_plate()
    mdl.create_explicit_step()
    mdl.create_assembly()
    mdl.constrain_b_plate()
    mdl.constrain_clamps()
    mdl.constrain_impactor()
    mdl.apply_symmetry()
    mdl.create_interactions_explicit()
    mdl.create_job()
    mdl.save_model()
