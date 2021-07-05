
from abaqusConstants import *
from sys import path
path.append('C:\\Python27\\Lib\\site-packages')
path.append('C:\\GitHub\\interlaced_model_creation\\mesoscale')
from interlaced_3d import Interlaced3D


class DCB(Interlaced3D):
    def __init__(self, model_name, dir_name=None):
        Interlaced3D.__init__(self, model_name, dir_name)

    def set_test_parameters(self, time, test_speed, output_intervals,
                            mesh_size):
        '''Set test parameters'''
        self.time = time
        self.test_speed = test_speed
        self.output_intervals = output_intervals
        self.mesh_size = mesh_size

    def create_beam(self, length, width, height, insert):
        self.length = length
        self.width = width
        self.height = height
        self.insert = insert
        sketch = self.model.ConstrainedSketch(
            name='Beam Profile', sheetSize=200.0)
        sketch.rectangle(point1=(0.0, 0.0), point2=(self.length, self.height))
        part = self.model.Part(
            name='Beam', dimensionality=THREE_D, type=DEFORMABLE_BODY)
        part.BaseSolidExtrude(sketch=sketch, depth=self.width)
        dPlane_id = part.DatumPlaneByPrincipalPlane(
            principalPlane=YZPLANE, offset=self.insert).id
        dPlane = part.datums[dPlane_id]
        part.PartitionCellByDatumPlane(datumPlane=dPlane, cells=part.cells[:])
        region = part.Set(cells=part.cells[:], name='Beam Region')
        part.SectionAssignment(
            region=region, sectionName='Tape-Elastic', offset=0.0,
            offsetType=MIDDLE_SURFACE, offsetField='',
            thicknessAssignment=FROM_SECTION)
        part.MaterialOrientation(
            region=region, orientationType=SYSTEM, stackDirection=STACK_3,
            angle=0.0, additionalRotationType=ROTATION_ANGLE,
            axis=AXIS_3, fieldName='', additionalRotationField='',
            localCsys=None)
        self.beam_part = part
        self.mesh_part_3d(self.beam_part, self.mesh_size)  # remesh

    def create_assembly(self):
        self.top_instance = self.assembly.Instance(
            name='Top', part=self.beam_part, dependent=ON)
        self.bottom_instance = self.assembly.Instance(
            name='Bottom', part=self.beam_part, dependent=ON)
        self.assembly.translate(
            instanceList=('Bottom', ), vector=(0.0, -self.height, 0.0))

    def apply_constraints(self):
        '''
        This method applies constraints to the assembly.
        '''
        # identify edges at top and bottom for coupling
        top_edge = self.top_instance.edges.findAt(
            ((0.0, self.height, self.width / 2.0), ))
        top_edge_surface = self.assembly.Set(
            edges=top_edge, name='Top Edge')
        bottom_edge = self.bottom_instance.edges.findAt(
            ((0.0, -self.height, self.width / 2.0), ))
        bottom_edge_surface = self.assembly.Set(
            edges=bottom_edge, name='Bottom Edge')

        # create reference points for coupling
        rf_point_1_id = self.assembly.ReferencePoint(
            point=(0.0, self.height + 1.0, self.width / 2.0)).id
        rf_point_1 = self.assembly.referencePoints[rf_point_1_id]
        rf_point_1_region = self.assembly.Set(
            referencePoints=(rf_point_1,), name='Top Reference Point')
        rf_point_2_id = self.assembly.ReferencePoint(
            point=(0.0, -self.height - 1.0, self.width / 2.0)).id
        rf_point_2 = self.assembly.referencePoints[rf_point_2_id]
        rf_point_2_region = self.assembly.Set(
            referencePoints=(rf_point_2,), name='Bottom Reference Point')

        self.model.Coupling(
            name='Top Coupling', controlPoint=rf_point_1_region,
            surface=top_edge_surface, influenceRadius=WHOLE_SURFACE,
            couplingType=KINEMATIC, localCsys=None, u1=ON, u2=ON, u3=ON,
            ur1=ON, ur2=ON, ur3=ON)
        self.model.Coupling(
            name='Bottom Coupling', controlPoint=rf_point_2_region,
            surface=bottom_edge_surface, influenceRadius=WHOLE_SURFACE,
            couplingType=KINEMATIC, localCsys=None, u1=ON, u2=ON,
            u3=ON, ur1=ON, ur2=ON, ur3=ON)

        self.model.SmoothStepAmplitude(
            name='Smoothing Amplitude', timeSpan=STEP,
            data=((0.0, 0.0), (1e-05, 1.0)))
        self.model.VelocityBC(
            name='Top BC', createStepName='Loading Step',
            region=rf_point_1_region, v1=0.0, v2=self.test_speed, v3=0.0,
            vr1=0.0, vr2=0.0, vr3=UNSET, amplitude='Smoothing Amplitude',
            localCsys=None, distributionType=UNIFORM, fieldName='')
        self.model.DisplacementBC(
            name='Bottom BC', createStepName='Loading Step',
            region=rf_point_2_region, u1=0.0, u2=0.0, u3=0.0, ur1=0.0, ur2=0.0,
            ur3=0.0, amplitude=UNSET, fixed=OFF, distributionType=UNIFORM,
            fieldName='', localCsys=None)

    def create_interactions(self):
        # create explicit general contact definition
        self.model.ContactExp(name='GC', createStepName='Initial')
        general_contact = self.model.interactions['GC']
        general_contact.includedPairs.setValuesInStep(
            stepName='Initial', useAllstar=ON)
        # define 'Tangential' as default behaviour
        general_contact.contactPropertyAssignments.appendInStep(
            stepName='Initial', assignments=((GLOBAL, SELF, 'Tangential'), ))
        # assign cohesive behaviour to contacting tape surfaces
        top_cohesive_face = self.top_instance.faces.findAt(
            ((self.insert + 1.0, 0.0, self.width / 2.0), ))
        bottom_cohesive_face = self.bottom_instance.faces.findAt(
            ((self.insert + 1.0, 0.0, self.width / 2.0), ))
        top_surface = self.assembly.Surface(
            side1Faces=top_cohesive_face, name='Top Cohesive')
        bottom_surface = self.assembly.Surface(
            side1Faces=bottom_cohesive_face, name='Bottom Cohesive')
        general_contact.contactPropertyAssignments.appendInStep(
            stepName='Initial',
            assignments=((top_surface, bottom_surface, 'Cohesive'), ))

    def create_materials(self, E11, E22, E33, nu12, nu13, nu23, G12, G13,
                         G23, density):
        mat = self.model.Material(name='Tape-Elastic')
        mat.Elastic(type=ENGINEERING_CONSTANTS,
                    table=((E11, E22, E33, nu12, nu13, nu23, G12, G13, G23), ))
        mat.Density(table=((density, ), ))
        self.model.HomogeneousSolidSection(
            name='Tape-Elastic', material='Tape-Elastic')

if __name__ == '__main__':

    # Material input data
    E11 = 116.6
    E22 = E33 = 7.231
    nu12 = nu13 = 0.339
    nu23 = 0.374
    G12 = G13 = 3.268
    G23 = 2.632
    density = 1.59e-6

    # Create model
    mdl = DCB('DCB-Model')
    mdl.set_test_parameters(5.0, 1.0, 100, 0.5)
    mdl.create_materials(E11, E22, E33, nu12, nu13, nu23, G12, G13,
                         G23, density)
    mdl.create_beam(125.0, 25.0, 3.0, 30.0)
    mdl.create_assembly()
    mdl.create_explicit_step()
    mdl.t_thickness = 3.0
    mdl.create_interaction_properties(mdl.mesh_size)
    mdl.create_interactions()
    mdl.apply_constraints()
    mdl.create_job()
    mdl.save_model()
