from abaqus import *
from abaqusConstants import *
from sys import path
path.append('C:\\Python27\\Lib\\site-packages')
path.append('C:\\GitHub\\interlaced_model_creation\\mesoscale')
from shapely.geometry import Polygon
from interlaced_3d_2 import Interlaced3D


class TensileModel(Interlaced3D):
    def __init__(self, model_name, dir_name=None):
        Interlaced3D.__init__(self, model_name, dir_name)

    def set_test_parameters(self, time, test_speed, output_intervals,
                            mesh_size, damage_mode, symmetry_mode):
        '''Set test parameters'''
        self.time = time  # simulation time in ms
        self.test_speed = test_speed  # crosshead velocity in mm/s
        self.output_intervals = output_intervals  # number of output intervals
        self.mesh_size = mesh_size  # mesh size in mm
        self.damage_mode = damage_mode  # damage mode: 'ON' or 'OFF'
        self.symmetry_mode = symmetry_mode  # 'BC' or 'GEOM'

    def apply_constraints(self):
        '''
        This method applies constraints to the assembly.

        Args:
            None

        Returns:
            None
        '''
        tol = 0.001
        if self.symmetry_mode == 'BC':
            z_min = 0.0
            z_max = self.l_thickness / 2.0
        else:
            z_min = -self.l_thickness / 2.0
            z_max = self.l_thickness / 2.0

        # identify faces at top and bottom for coupling
        top_faces_box = (
            self.x_min - tol, self.y_max - tol, z_min - tol,
            self.x_max + tol, self.y_max + tol, z_max + tol)
        top_faces = self.get_faces_by_bounding_box(top_faces_box)
        top_faces_surface = self.assembly.Surface(
            side1Faces=top_faces, name='Top Faces')
        bottom_faces_box = (
            self.x_min - tol, self.y_min - tol, z_min - tol,
            self.x_max + tol, self.y_min + tol, z_max + tol)
        bottom_faces = self.get_faces_by_bounding_box(bottom_faces_box)
        bottom_faces_surface = self.assembly.Surface(
            side1Faces=bottom_faces, name='Bottom Faces')

        # create reference points for coupling
        rf_point_1_id = self.assembly.ReferencePoint(
            point=(0.0, y_max + 5.0, 0.0)).id
        rf_point_1 = self.assembly.referencePoints[rf_point_1_id]
        rf_point_1_region = self.assembly.Set(
            referencePoints=(rf_point_1,), name='Reference Point 1')
        rf_point_2_id = self.assembly.ReferencePoint(
            point=(0.0, y_min - 5.0, 0.0)).id
        rf_point_2 = self.assembly.referencePoints[rf_point_2_id]
        rf_point_2_region = self.assembly.Set(
            referencePoints=(rf_point_2,), name='Reference Point 2')

        self.model.Coupling(
            name='Top Coupling', controlPoint=rf_point_1_region,
            surface=top_faces_surface, influenceRadius=WHOLE_SURFACE,
            couplingType=KINEMATIC, localCsys=None, u1=ON, u2=ON, u3=ON,
            ur1=ON, ur2=ON, ur3=ON)

        self.model.Coupling(
            name='Bottom Coupling', controlPoint=rf_point_2_region,
            surface=bottom_faces_surface, influenceRadius=WHOLE_SURFACE,
            couplingType=KINEMATIC, localCsys=None, u1=ON, u2=ON,
            u3=ON, ur1=ON, ur2=ON, ur3=ON)

        self.model.SmoothStepAmplitude(
            name='Smoothing Amplitude', timeSpan=STEP,
            data=((0.0, 0.0), (1e-05, 1.0)))

        self.model.DisplacementBC(
            name='Bottom Surface BC', createStepName='Loading Step',
            region=rf_point_2_region, u1=0.0, u2=0.0, u3=0.0, ur1=0.0, ur2=0.0,
            ur3=0.0, amplitude=UNSET, fixed=OFF, distributionType=UNIFORM,
            fieldName='', localCsys=None)

        if self.damage_mode == 'ON':    
            self.model.VelocityBC(
                name='Top Surface BC', createStepName='Loading Step',
                region=rf_point_1_region, v1=0.0, v2=self.test_speed, v3=0.0,
                vr1=0.0, vr2=0.0, vr3=0.0, amplitude='Smoothing Amplitude',
                localCsys=None, distributionType=UNIFORM, fieldName='')
        else:
            displacement = self.test_speed * self.time
            self.model.DisplacementBC(
                name='Top Surface BC', createStepName='Loading Step',
                region=rf_point_1_region, u1=0.0, u2=displacement, u3=0.0,
                ur1=0.0, ur2=0.0, ur3=0.0, amplitude=UNSET, fixed=OFF,
                distributionType=UNIFORM, fieldName='', localCsys=None)

        if self.symmetry_mode == 'BC':
            symmetric_faces_box = (
                self.x_min - tol, self.y_min - tol, -tol,
                self.x_max + tol, self.y_max + tol, tol)
            symmetric_faces = self.get_faces_by_bounding_box(
                symmetric_faces_box)
            symmetric_faces_set = self.assembly.Set(
                faces=symmetric_faces, name='Symmetry Faces')
            self.model.ZsymmBC(
                name='Symmetry', createStepName='Loading Step',
                region=symmetric_faces_set, localCsys=None)


if __name__ == '__main__':

    # Simulation parameters
    time = 5.0  # duration to simulate [ms]
    output_intervals = 50  # requested field output intervals
    crosshead_velocity = 0.5
    damage_mode = 'OFF'
    symmetry_mode = 'GEOM'

    # Material parameters
    tape_angles = (0, 45, 90, -45)  # define angles of tapes
    tape_widths = 10.0
    tape_spacing = 3  # number of gaps between tapes in interlacing pattern
    cured_ply_thickness = 0.18  # cured ply thickness e.g. 0.18125
    undulation_ratio = 0.09  # ratio of undulation amplitude to length
    number_of_plies = 8  # symmetric 4 ply laminate

    # Mesh sizes
    mesh_size = 1.0

    # RVE dimensions
    x_min = -20.0
    x_max = 20.0
    y_min = -50.0
    y_max = 50.0
    specimen_size = Polygon([(x_min, y_min), (x_max, y_min), (x_max, y_max),
                             (x_min, y_max)])

    # Create model
    mdl = TensileModel('QI_implicit_3')
    mdl.create_implicit_step()
    mdl.set_test_parameters(time, crosshead_velocity, output_intervals,
                            mesh_size, damage_mode, symmetry_mode)
    mdl.set_specimen_parameters(tape_angles, tape_widths, tape_spacing,
                                cured_ply_thickness, undulation_ratio,
                                number_of_plies)
    mdl.create_materials()
    paths_f, angles_f, grid_f = mdl.create_tape_paths(specimen_size)
    mdl.create_parts(grid_f)
    mdl.create_sequences()
    mdl.create_cohesive_interactions(tie=True)
    mdl.apply_constraints()
    mdl.create_job()
    mdl.save_model()
